function [X] = langles_to_vecs3d(m, fids, angles)
    % m - input mesh
    % fids - list of face ids
    % angles - angles corresponding to given faces. Each angle is 
    % relative to the local frame of that face.
    %
    % Output:
    %   X - the output 3d coordinates. (length(fids) x 3 matrix)
    [local_frames, ~] = create_local_frames(m);
    nf = m.nF;
    X = zeros(length(fids), 3);
    for i = 1:length(fids)
        fid = fids(i);
        theta = angles(i);
        frame = [local_frames(fid, :); ...
                 local_frames(fid+nf, :)];
        g_vec = [cos(theta), sin(theta)] * frame;
        X(i, :) = g_vec;
    end
end

