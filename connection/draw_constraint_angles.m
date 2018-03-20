function draw_constraint_angles(m, constrained_faces, constraints)
    % function draw_constraints(m, constrained_faces, constraints)
    %
    % Input:
    %   m - mesh object
    %   constrained_faces - nC x 1 ids of the constrained faces
    %   constraint_vectors - nC x 3 matrix of constraints
    
    %if size(constraints, 2) ~= 3
        [local_frames, ~] = create_local_frames(m);
        nC = length(constrained_faces);
        constraint_vectors = zeros(nC, 3);
        for k = 1:nC
            fid = constrained_faces(k);
            t = constraints(k);
            nF = m.nF;
            frame = [local_frames(fid, :); ...
                     local_frames(fid+nF, :)];
            vec = [cos(t), sin(t)] * frame;
            constraint_vectors(k, :) = vec;
        end
    %else
    %    constraint_vectors = constraints;
    %end

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

