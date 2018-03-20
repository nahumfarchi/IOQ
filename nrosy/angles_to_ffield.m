function [ffield] = angles_to_ffield(thetas, local_frames, N, fids_list)
    % function [NROSY] = angles_to_nrosy(thetas, local_frames, N)
    %
    % Given angles relative to the given frame of reference,
    % output an nrosy field in global coordinates.
    %
    % Input:
    %   thetas - nFx1 vector of angles relative to local_frames
    %   local_frames - 2*nFx3 matrix of local frames, first e1, then e2.
    %   N - nrosy field degree
    %
    % Output:
    %   NROSY - N*nFx3 matrix of tangent vectors. First nF rows are the
    %   first tangent vector of each face, second nF rows are the second
    %   tangent vector, etc.
    
    nF = size(local_frames, 1) / 2;
    if nargin < 4
        fids_list = 1:nF;
    end
    ffield = zeros(N*length(fids_list), 3);
    %for fid = fids_list
    for k = 1:length(fids_list)
        fid = fids_list(k);
        % assume frame is orthogonal
        frame = [local_frames(fid, :); ...
                 local_frames(fid+nF, :)];
        for i = 0:(N-1)
            t = thetas(k) + i*2*pi/N;
            vec = [cos(t), sin(t)] * frame;
            %vec = frame * [cos(t); sin(t); 0];
            ffield(k+i*length(fids_list), :) = vec;
        end
    end

end

