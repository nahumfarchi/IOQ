function DF = local_angles_to_gvectors(m, thetas, local_frames, N)
    % function DF = local_angles_to_gvectors(m, thetas, local_frames, N)
    %
    % Input:
    %   mesh object
    %   theta - F x 1 vector of angles
    %   local_frames - 2*F x 3 matrix of vector frames. First F rows are
    %   for e1, second F rows are for e2.
    %   
    % Output:
    %   DF - (N*F) x 3 direction field matrix. First F rows are v1, second
    %   F rows are v2, etc.
    DF = zeros(N*m.nF, 3);
    for fid = 1:m.nF
        if isnan(thetas(fid))
            continue;
        end
        frame = [local_frames(fid, :); ...
                 local_frames(fid+m.nF, :)];
        for i = 0:(N-1)
            t = thetas(fid) + i*2*pi/N;
            vec = [cos(t), sin(t)] * frame;
            DF(fid+i*m.nF, :) = vec;
        end
    end
    DF = DF(~isnan(thetas), :);
end