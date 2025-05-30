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
            addOptional(p,'Dock',1);
            parse(p,varargin{:});
            
            patch('faces',M.triangles,'vertices',M.vertices, ...
                  'FaceColor',p.Results.FaceColor, ...
                  'EdgeColor',p.Results.EdgeColor, 'FaceAlpha',0.9);
              
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
            addOptional(p,'Dock',1);
            addOptional(p, 'Colormap', 'jet');
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

            patch('faces',M.triangles,'vertices',M.vertices, ...
                  CW,f, ...
                  'FaceColor',FC, ...
                  'EdgeColor',p.Results.EdgeColor); 
            
            cameratoolbar; cameratoolbar('SetCoordSys','none');
            axis equal; axis off;
            colorbar; caxis(p.Results.Caxis); 
            colormap(p.Results.Colormap);
%             view(p.Results.View);
            if p.Results.Dock == 1; set(gcf, 'WindowStyle', 'docked'); end;
            
%             view([0 1 0]);
        end
        
        function vf(M,vf,varargin)
            vf = reshape(vf,[],3);
            
            p = inputParser;
            addOptional(p,'FaceColor','w');
            addOptional(p,'EdgeColor','k');
            addOptional(p,'f',MESH.normv(vf));
            addOptional(p,'Color','k');
            addOptional(p,'View',[0 1 0]);
            addOptional(p,'Dock',0);
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
            elseif length(vf) == M.nv
                Xm = M.vertices;
            end
            
            a = .05; vf = a*vf;
            quiver3(Xm(:,1),Xm(:,2),Xm(:,3), ...
                    vf(:,1),vf(:,2),vf(:,3), ...
                    'Color',p.Results.Color);
                
            cameratoolbar; cameratoolbar('SetCoordSys','none');
            axis equal; axis off; 
            colormap(jet);
            if p.Results.Dock == 1; set(gcf, 'WindowStyle', 'docked'); end;
            
%             view([0 1 0]);
            
            hold off;
        end
    end
    
end

