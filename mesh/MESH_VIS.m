classdef MESH_VIS
    %MESH_VIS Summary of this class goes here
    %   Detailed explanation goes here
    % TODO: mesh,vf
        
    properties
    end
    
    methods (Static)
        
        function mesh(M,varargin)
            p = inputParser;
            addOptional(p,'FaceColor','w');
            addOptional(p,'EdgeColor','k');
            addOptional(p,'CW',[]);
            addOptional(p,'UVt',[]);
            addOptional(p,'FUVt',[]);
            addOptional(p,'Texture','cross.png');
            addOptional(p,'SingularPts',0);
            addOptional(p,'SingularPtsR',20);
            parse(p,varargin{:});
            
            if isempty( p.Results.UVt )
                if any(strcmp(properties(M), 'triangles')) == 1
                    patch('faces',M.triangles,'vertices',M.vertices, ...
                          'FaceColor',p.Results.FaceColor, ...
                          'EdgeColor',p.Results.EdgeColor);
                else
                    if isempty(p.Results.CW) == 1
                        patch('faces',M.quads,'vertices',M.vertices, ...
                              'FaceColor',p.Results.FaceColor, ...
                              'EdgeColor',p.Results.EdgeColor);
                    else
                        
                        CW = 'CData';
                        if size(p.Results.CW,2) == 3; CW = 'FaceVertexCData'; end;
                        
                        patch('faces',M.quads,'vertices',M.vertices, ...
                              CW,p.Results.CW,...
                              'FaceColor','interp', ...
                              'EdgeColor',p.Results.EdgeColor,...
                              'LineWidth',2);
                    end
                    
                    if p.Results.SingularPts == 1 
                        MESH_VIS.singular_pts( M, p.Results.SingularPtsR );
%                         camlight; lighting phong; material dull;
                    end
                end
            else
                I = imread(p.Results.Texture);
                
                patcht2(M.triangles,M.vertices,...
                        p.Results.FUVt,p.Results.UVt,I);
            end
              
            cameratoolbar; cameratoolbar('SetCoordSys','none'); 
            axis equal; axis off;
        end
        
        function [ M2 ] = nc_mesh(M)
            M2.vertices = ( M.vertices(M.edges(:,1),:) + ...
                            M.vertices(M.edges(:,2),:) )/2;
            M2.triangles = abs(M.t2e);

            M2.nf = M.nf;
            M2.nv = M.ne;

            M2.om = M;
            M2.om.vertices = M.vertices-M.Nv*M.mel/10;
        end
        
        function func(M,f,varargin)
            p = inputParser;
            addOptional(p,'EdgeColor','none');
            addOptional(p,'Caxis','auto');
            addOptional(p,'View',[0 1 0]);
            addOptional(p,'Dock',0);
            addOptional(p,'Colormap','jet');
            addOptional(p,'OpenGL',0);
            addOptional(p,'LightPos',[0 0 1]);
            parse(p,varargin{:});
            
            szf = size(f);
            if szf(1) == szf(2) && norm(f - diag(diag(f)),'fro') < 1e-8
                f = diag(f);
            end
            
            FC = 'w';
            if szf(1) == M.nv; FC = 'interp';
            elseif szf(1) == M.nf; FC = 'flat'; 
            else
                M = MESH_VIS.nc_mesh( M );
                patch('faces',M.om.triangles,'vertices',M.vertices, ...
                      'FaceColor',FC, ...
                      'EdgeColor',p.Results.EdgeColor);
                FC = 'interp';
            end;
            
            CW = 'CData';
            if szf(2) == 3; CW = 'FaceVertexCData'; end;
            
            if any(strcmp(properties(M), 'triangles')) == 1
                patch('faces',M.triangles,'vertices',M.vertices, ...
                      CW,f, ...
                    'FaceColor',FC, ...
                    'EdgeColor',p.Results.EdgeColor); 
            else
                patch('faces',M.quads,'vertices',M.vertices, ...
                      CW,f, ...
                    'FaceColor',FC, ...
                    'EdgeColor',p.Results.EdgeColor); 
            end
            
            cameratoolbar; cameratoolbar('SetCoordSys','none');
            axis equal; axis off;
            colorbar; caxis(p.Results.Caxis); 
            colormap(p.Results.Colormap);

            if p.Results.OpenGL == 1
                camlight('headlight'); 
                lighting phong;
                material dull;
%                 light('Position',p.Results.LightPos,'Style','infinite');
            end
                
            if p.Results.Dock == 1; set(gcf, 'WindowStyle', 'docked'); end;
        end
        
        function vf(M,vf,varargin)
            vf = reshape(vf,[],3);
            
            p = inputParser;
            addOptional(p,'FaceColor','w');
            addOptional(p,'EdgeColor','k');
            addOptional(p,'f',MESH.normv(vf));
            addOptional(p,'Color','k');
            addOptional(p,'Dock',0);
            addOptional(p,'nRosy',1);
            parse(p,varargin{:});
            
            hold on;
            if isempty(p.Results.f) == 1            
                patch('faces',M.triangles,'vertices',M.vertices, ...
                      'FaceColor',p.Results.FaceColor, ...
                      'EdgeColor',p.Results.EdgeColor);
            else
                MESH_VIS.func(M,p.Results.f);
            end
            
            if length(vf) == M.nf
                X = M.vertices; T = M.triangles;
                Xm = (X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:))/3;
            elseif length(vf) == M.ne
                X = M.vertices; E = M.edges;
                Xm = X(E(:,1),:)+X(E(:,2),:)/2;
            end
            
            scale = .5;
            R = speye(3*M.nf,3*M.nf);
            for i = 1:(4/p.Results.nRosy); R = M.R*R; end;

            % works only for n=2/4
            for i = 1:p.Results.nRosy
                quiver3( Xm(:,1), Xm(:,2), Xm(:,3), ...
                         vf(:,1), vf(:,2), vf(:,3), ...
                         scale, 'Color', p.Results.Color );

                 vf = reshape(R*vf(:),[],3);
            end
                
            cameratoolbar; cameratoolbar('SetCoordSys','none');
            axis equal; axis off; 
            colormap(jet);
            if p.Results.Dock == 1; set(gcf, 'WindowStyle', 'docked'); end;
            
%             view([0 1 0]);
            
            hold off;
        end
        
        function set_camera(ca,cam)
            set(ca, 'PlotBoxAspectRatio',cam.pba);
            set(ca, 'DataAspectRatio',cam.dar);
            set(ca, 'CameraViewAngle',cam.cva);
            set(ca, 'CameraUpVector',cam.cuv);
            set(ca, 'CameraTarget',cam.ct);
            set(ca, 'CameraPosition',cam.cp);
        end
        
        function cam = get_camera(ca)
            cam.pba = get(ca, 'PlotBoxAspectRatio');
            cam.dar = get(ca, 'DataAspectRatio');
            cam.cva = get(ca, 'CameraViewAngle');
            cam.cuv = get(ca, 'CameraUpVector');
            cam.ct = get(ca, 'CameraTarget');
            cam.cp = get(ca, 'CameraPosition');
        end
        
        function singular_pts( M, r )
            X = M.vertices; 
            
            SI3 = find(M.valence == 3);
            SI5 = find(M.valence == 5);
            
            [SX,SY,SZ] = sphere; 
            a = max(max(X)-min(X))*r;
            SX = SX/a; SY = SY/a; SZ = SZ/a;

            hold on;
            for i = 1:numel(SI3)
                xs = X(SI3(i),:);
                surf(SX+xs(1),SY+xs(2),SZ+xs(3),'FaceColor','r',...
                    'EdgeColor','none','FaceAlpha',1);
            end
            
            for i = 1:numel(SI5)
                xs = X(SI5(i),:);
                surf(SX+xs(1),SY+xs(2),SZ+xs(3),'FaceColor','b',...
                    'EdgeColor','none','FaceAlpha',1);
            end
            hold off;
        end
        
        function cm2 = interp_cm( cm, sz )
            n = numel(cm(:,1));
            xq = linspace(0,1,sz);
            cm2 = zeros(sz,3);
            for i = 1:3
                cm2(:,i) = interp1(cm(:,1),cm(:,i+1),xq);
            end
        end
    end
    
end

