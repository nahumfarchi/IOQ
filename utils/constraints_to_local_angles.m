function [angles] = constraints_to_local_angles(m, local_frames, constrained_faces, constraint_vectors)
    % Given constraints vectors in global coordinates, convert to local
    % angles. For each constraint, project the constraint vector unto the local frame of the
    % corresponding face and return the angle.
    %
    % Input:
    %   local_frames - (2*F) x 3 matrix of local frames per face.
    %   constrained_faces - nc x 1 vector with the ids of the constrained
    %   faces
    %   constraint_vectors - nc x 3 matrix of constraint vectors in global
    %   coordinates
    %
    % Output:
    %   angles - the constraints angles relative to the corresponding local
    %   frames.

    nc = length(constrained_faces); % number of constraints
    angles = zeros(nc, 1);
    for i = 1:length(constrained_faces)
        fid = constrained_faces(i);
        v = constraint_vectors(i, :);
        frame = [local_frames(fid, :);
                 local_frames(fid+m.nF, :)];
        v_proj = frame * v';
        angles(i) = atan2(v_proj(2), v_proj(1));
    end
   
end

