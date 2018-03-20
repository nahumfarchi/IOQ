function [] = draw_mesh(m, color, face_alpha, edge_alpha)
% function [] = draw_mesh(m, color, face_alpha, edge_alpha)
%
% Draw the given mesh.
%
% Input:
%   m - mesh object
%   color - color string or 1x3 vector
%   face_alpha, edge_alpha (optional)

    if nargin < 3
        face_alpha = 1;
    end
    if nargin < 4
        edge_alpha = 0.1;
    end

    patch('Faces', m.F, 'Vertices', m.V, 'FaceColor', color, 'FaceAlpha', face_alpha, 'EdgeAlpha', edge_alpha)

    view(3);
    ax = gca;
    ax.Clipping = 'off';

    camproj('perspective');
    axis square; 
    %axis off;

    axis tight;
    axis equal;
    cameramenu;

end

