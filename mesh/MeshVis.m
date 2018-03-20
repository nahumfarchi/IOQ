classdef MeshVis
    %MESH_VIS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods (Static)
        
        function plot(m, varargin)
            p = inputParser;
            addOptional(p, 'func', []);
            addOptional(p, 'FaceAlpha', 1.0);
            addOptional(p, 'EdgeAlpha', 0.2);
            addOptional(p, 'nrosy', []);
            addOptional(p, 'nrosy_colors', {'b', 'b', 'b', 'b'});
            addOptional(p, 'scale', 1);
            addOptional(p, 'colormap', 'jet');
            addOptional(p, 'ConstrainedFaces', []);
            addOptional(p, 'ConstraintAngles', []);
            addOptional(p, 'ConstraintVectors', []);
            addOptional(p, 'LabelEdges', false);
            addOptional(p, 'LabelFaces', false);
            addOptional(p, 'LabelVertices', false);
            addOptional(p, 'LocalFrames', false);
            addOptional(p, 'Title', []);
            parse(p, varargin{:});
            
            func = p.Results.func;
            if isempty(func)
                func = ones(m.nF, 1);
                FaceColor = 'white';
            elseif length(func) == m.nV
                FaceColor = 'interp';
            elseif length(func) == m.nF
                FaceColor = 'flat';
            else
                error('Bad func size')
            end
            func = func(:);
            
            constrained_faces = p.Results.ConstrainedFaces;
            constraint_angles = p.Results.ConstraintAngles;
            constraint_vectors = p.Results.ConstraintVectors;
            
            nrosy = p.Results.nrosy;
            hold on
            if ~isempty(nrosy)
                % Plot singularities
                if isfield(nrosy, 'S')
                    S = nrosy.S;
                    for j = 1:size(S,1)
                        fid = S(j, 1);
                        if S(j, 2) > 0
                            H = plot3( m.V(fid,1), m.V(fid,2), m.V(fid,3), 'r.' );
                        else
                            H = plot3( m.V(fid,1), m.V(fid,2), m.V(fid,3), 'g.' );
                        end
                        set( H, 'MarkerSize', 40 );
                    end
                end
                
                % Plot nrosy field
                ffield = nrosy.ffield;
                degree = nrosy.degree;
                nrosy_colors = p.Results.nrosy_colors;
                scale = p.Results.scale;
                for j = 0:degree-1
                    inds = 1+j*m.nF : m.nF+j*m.nF;
                    m.drawFaceField(ffield(inds, :), ...
                        'color', nrosy_colors{j+1}, ...
                        'AutoScale', 'on', ...
                        'AutoScaleFactor', scale)
                end
            end
                     
            patch('faces', m.F, ...
                  'vertices', m.V, ...
                  'FaceVertexCData', func, ...
                  'FaceColor', FaceColor, ...
                  'FaceAlpha', p.Results.FaceAlpha, ...
                  'EdgeAlpha', p.Results.EdgeAlpha);
              
            if ~isempty(constraint_angles)
                draw_constraint_angles(m, constrained_faces, constraint_angles); 
            elseif ~isempty(constraint_vectors)
                P = (m.V(m.F(constrained_faces, 1), :) + ...
                m.V(m.F(constrained_faces, 2), :) + ...
                m.V(m.F(constrained_faces, 3), :)) / 3;
                x = P(:, 1);
                y = P(:, 2);
                z = P(:, 3);
                u = constraint_vectors(:, 1);
                v = constraint_vectors(:, 2);
                w = constraint_vectors(:, 3);
                quiver3(x, y, z, u, v, w, 'Marker', 'o', 'AutoScale', 'on', 'AutoScaleFactor', m.avg_length, 'color', 'g', 'LineWidth', 2);
            end
            
            if numel(p.Results.LabelFaces) == m.nF
                m.labelFaces(p.Results.LabelFaces, '%.2g');
            elseif p.Results.LabelFaces
                m.labelFaces();
            end
            if numel(p.Results.LabelEdges) == m.nE
                m.labelEdges(p.Results.LabelEdges, '%.2g');
            elseif p.Results.LabelEdges
                m.labelEdges();
            end
            if numel(p.Results.LabelVertices) == m.nV
                m.labelVertices(p.Results.LabelVertices, '%.2g');
            elseif p.Results.LabelVertices
                m.labelVertices();
            end
            
            if p.Results.LocalFrames
                [local_frames, ~] = create_local_frames(m);
                draw_local_frames(m, local_frames, p.Results.scale)
            end
              
            hold off
            
            if ~isempty(p.Results.Title)
                title(p.Results.Title, ...
                        'FontSize',25,...
                        'interpreter','latex');   
            end
            
            colormap(p.Results.colormap);
              
            view(3);
            ax = gca;
            ax.Clipping = 'off';
            
            camproj('perspective');
            axis square; 
            axis off;

            axis tight;
            axis equal;
            cameramenu;
        end
        
        function wfigs(filename, meshes, varargin)
            % Usage: TODO
            %   
            
            p = inputParser;
            addOptional(p,'Titles',[]);
            addOptional(p,'Colormap','jet');
            addOptional(p,'CAxis','auto');
            addOptional(p,'CAxisV',[]);
%             addOptional(p,'Plot','func');
            addOptional(p,'Montage',1);
            addOptional(p,'View',[]);
            addOptional(p,'Zoom',1);
            addOptional(p,'Resolution',1024);
            addOptional(p,'Camera',[]);
            addOptional(p,'OpenGL',0);
            %addOptional(p,'LightPos',[0 0 1]);
%             addOptional(p,'UVt',[]);
%             addOptional(p,'FUVt',[]);
%             addOptional(p,'SingularPts',0);
%             addOptional(p,'SingularPtsR',20);
            addOptional(p,'AspectRatio',1);
%             addOptional(p,'EdgeColor','k');
            addOptional(p, 'Funcs', {});
            addOptional(p, 'FaceAlpha', 1.0);
            addOptional(p, 'EdgeAlpha', 0.2);
            addOptional(p, 'Nrosy', {});
            addOptional(p, 'NrosyColors', {'b', 'b', 'b', 'b'});
            addOptional(p, 'Scale', 1);
            addOptional(p, 'CloseFigs', false);
            addOptional(p, 'ConstrainedFaces', []);
            addOptional(p, 'ConstraintAngles', []);
            addOptional(p, 'ConstraintVectors', []);
            addOptional(p, 'OutFolder', '.');
            addOptional(p, 'Save', true);
            addOptional(p, 'LabelEdges', false);
            addOptional(p, 'LabelFaces', false);
            addOptional(p, 'LabelVertices', false);
            addOptional(p, 'LocalFrames', false);
            addOptional(p, 'Rows', 2);
            parse(p,varargin{:});
            
            funcs = p.Results.Funcs;
            func = [];
            nrosys = p.Results.Nrosy;
            nrosy = [];
            if isempty(funcs) && isempty(nrosys)
               nFigs = 1;
               func = {};
            elseif ~isempty(funcs)
               nFigs = numel(funcs);
            elseif ~isempty(nrosys)
               nFigs = numel(nrosys);
            end
            if ~isempty(funcs) && ~isempty(nrosys)
               assert(numel(funcs) == numel(nrosys))
            end
            close_figs = p.Results.CloseFigs;
            view_angle = p.Results.View;
            
            out_folder = p.Results.OutFolder;
            [~, name_only, ~] = fileparts(filename);
            filename = fullfile(out_folder, name_only);    
            
            for i = 1:nFigs
                hdl = figure; 
                
                if ~isempty(funcs)
                    func = funcs{i};
                    colorbar
                end
                if ~isempty(nrosys)
                    nrosy = nrosys{i};
                end
                
                if numel(meshes) == nFigs && nFigs > 1
                    mesh = meshes{i};
                else
                    mesh = meshes;
                end
                if numel(p.Results.LabelFaces) == nFigs && nFigs > 1
                    label_faces = p.Results.LabelFaces{i};
                else
                    label_faces = p.Results.LabelFaces;
                end
                if numel(p.Results.LabelEdges) == nFigs && nFigs > 1
                    label_edges = p.Results.LabelEdges{i};
                else
                    label_edges = p.Results.LabelEdges;
                end
                if numel(p.Results.LabelVertices) == nFigs && nFigs > 1
                    label_vertices = p.Results.LabelVertices{i};
                else
                    label_vertices = p.Results.LabelVertices;
                end
                if numel(p.Results.View) == nFigs && nFigs > 1
                    view_angle = p.Results.View{i};
                end
                
                % Some weirdness to make sure that the title is not cut off
                sp3 = subplot(1,1,1);
                if numel(p.Results.Titles) > 0
                    title(p.Results.Titles{i}, ...
                        'FontSize',25,...
                        'interpreter','latex');
                elseif isfield(nrosy, 'title')
                    title(nrosy.title, ...
                        'FontSize',25,...
                        'interpreter','latex');
                end
                drawnow % force calculating the position *after* inserting the title
                ph = sp3.Position; % get the desired position
                sp3.delete % remove the axes
                sp3 = subplot(1,1,1);
                
                MeshVis.plot(mesh, ...
                    'func', func, ...
                    'FaceAlpha', p.Results.FaceAlpha, ...
                    'EdgeAlpha', p.Results.EdgeAlpha, ...
                    'nrosy', nrosy, ...
                    'nrosy_colors', p.Results.NrosyColors, ...
                    'scale', p.Results.Scale, ...
                    'colormap', p.Results.Colormap, ...
                    'ConstrainedFaces', p.Results.ConstrainedFaces, ...
                    'ConstraintAngles', p.Results.ConstraintAngles, ...
                    'ConstraintVectors', p.Results.ConstraintVectors, ...
                    'LabelFaces', label_faces, ...
                    'LabelEdges', label_edges, ...
                    'LabelVertices', label_vertices, ...
                    'LocalFrames', p.Results.LocalFrames);
                
                if p.Results.OpenGL == 1
                    lighting phong;
                    material dull;
%                   light('Position',p.Results.LightPos,...
%                         'Style','local','Color',[.8,.8,.8]);
                    camlight('headlight');
                end
                
                if ~isempty(view_angle)
                    view(view_angle);
                end
                
                if numel(p.Results.Titles) > 0
                    title(p.Results.Titles{i}, ...
                        'FontSize',25,...
                        'interpreter','latex');
                elseif isfield(nrosy, 'title')
                    title(nrosy.title, ...
                        'FontSize',25,...
                        'interpreter','latex');
                end
                
                sp3.Position = ph; % set the axes to the right height.
                % End of weirdness
                
                if p.Results.Save == 1
                    pngname = sprintf('%s_%04d.png', filename, i);

                    dpi = get(0, 'ScreenPixelsPerInch');
                    in = p.Results.Resolution/dpi;
                    
                    fig = gcf;
                    fig.PaperUnits = 'inches';
                    fig.PaperPosition = [0 0 p.Results.AspectRatio*in in];
                    fig.PaperPositionMode = 'manual';
                                      
                    print(pngname,'-dpng','-r0')
                    if close_figs
                        close(hdl);
                    end
                end
            end
            
            % montage
            if p.Results.Montage == 1 && p.Results.Save == 1
                if numel(funcs) <= 1 && numel(nrosys) <= 1
                   rows = 1;
                   cols = 1;
                elseif numel(funcs) > 1
                   rows = p.Results.Rows;
                   cols = ceil(length(funcs) / rows);
                elseif numel(nrosys) > 1
                   rows = p.Results.Rows;
                   cols = ceil(length(nrosys) / rows);
                else
                   rows = 1;
                   cols = nFigs;
                end
                %cols = size(F,2) / rows;

                A = [];
                k = 1;
                for i = 1:rows
                    b = [];
                    for j = 1:cols
                        pngname = sprintf('%s_%04d.png',filename,k);
                        if ~exist(pngname, 'file')
                            a = 255*ones(whd);
                        else
                            a = imread(pngname);
                            whd = size(a);
                        end
                        b = cat(2,b,a);
                        k = k + 1;
                    end
                    A = cat(1,A,b);
                end
                imwrite(A,[filename '.png']);
                delete( sprintf('%s_*.png',filename) );
                if close_figs
                    close all;
                end
            end

        end
        
     end
    
end

