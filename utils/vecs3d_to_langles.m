function [langles] = vecs3d_to_langles(m, fids, vecs3d)
    % Given vectors in global coordinates, convert to local
    % angles.
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

    [local_frames, ~] = create_local_frames(m);
    nc = length(fids);
    langles = zeros(nc, 1);
    for i = 1:length(fids)
        fid = fids(i);
        v = vecs3d(i, :);
        frame = [local_frames(fid, :);
                 local_frames(fid+m.nF, :)];
        v_proj = frame * v';
        langles(i) = atan2(v_proj(2), v_proj(1));
    end
end