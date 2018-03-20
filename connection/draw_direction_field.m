function draw_direction_field(m, DF, I, scale, cf_colors, face_alpha, edge_alpha, periods, thetas)
    % function draw_direction_field(m, DF, I, scale, cf_colors, face_alpha, edge_alpha, periods, thetas)
    %
    % Draw the given direction field on the mesh.
    %
    % Input:
    %   mesh object
    %   DF - (N*F) x 3 matrix with direction field vectors. First F rows 
    %   are v1, second F rows are v2, etc
    %   I - singularity indices    
    %   scale - scales factor
    %   cf_colors - eg {'r', 'b', 'b', 'b'}
    %   periods - if given, will label the edges with the periods
    %   thetas - if given, will label the faces with the angles
    
    N = size(DF, 1) / m.nF;
    label_periods = true;
    label_faces = true;
    if nargin < 3, I = []; end
    if nargin < 4, scale = 2*m.avg_length; end
    if nargin < 5, cf_colors = {'b', 'b', 'b', 'b'}; end
    if nargin < 6, face_alpha = 1.0; end
    if nargin < 7, edge_alpha = 0.2; end
    if nargin < 8, label_periods = false; end
    if nargin < 9, label_faces = false; end
    
    %m.draw(ones(m.nF, 1), 'EdgeAlpha', 0.1);
    draw_mesh(m, 'white', face_alpha, edge_alpha)
    
    hold on;
    
    for j = 0:N-1
        %DF = self.field(1+i*m.nF:m.nF+i*m.nF, :);
        ind = 1+j*m.nF : m.nF+j*m.nF;
        m.drawFaceField(DF(ind, :), 'color', cf_colors{j+1}, 'AutoScale', 'on', 'AutoScaleFactor', scale)
    end
    
    % plot singularities
    for j = 1:size(I,1)
        fid = I(j, 1);
        H = plot3( m.V(fid,1), m.V(fid,2), m.V(fid,3), 'r.' );
        set( H, 'MarkerSize', 40 );
    end
    
    if label_periods
        for eid = 1:m.nE
            %vid1 = m.EVAdj(eid, 1);
            %vid2 = m.EVAdj(eid, 2);
            v1 = m.V(m.EVAdj(eid, 1), :);
            v2 = m.V(m.EVAdj(eid, 2), :);
            p = (v1 + v2) / 2;
            %p = (m.V(vid1,:) + m.V(vid2,:)) ./ 2;
            %p = p + 0.1 * (m.V(vid2,:) - m.V(vid1,:));
            
            periods(abs(periods) < 1e-8) = 0;
            %if periods(eid) > 1e-8 || periods(eid) < -1e-8
                text(p(1), p(2), p(3), num2str(periods(eid)), 'FontSize', 8) %['e', num2str(eid)])
            %end
        end
    end
    
    if label_faces
        for fid = 1:m.nF
            v1 = m.V(m.F(fid, 1), :);
            v2 = m.V(m.F(fid, 2), :);
            v3 = m.V(m.F(fid, 3), :);
            p = (v1 + v2 + v3) / 3 + 0.05*m.avg_length*(v1-v2)/norm(v1-v2);
            text(p(1), p(2), p(3), num2str(thetas(fid), '%.2f'), 'FontSize', 8, 'color', 'r')
        end
    end

    axis equal;
    axis off;
    view([-180,0]);
    zoom(1);
    rotate3d;
    hold off;
end