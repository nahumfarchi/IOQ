function [langle] = gangles_to_local_angles(theta, local_frames, N, fids_list)
    % function [gangle] = gangles_to_local_angles(theta, local_frames, N, fids_list)
    %
    % Given angles relative XY axis, return angles in local frames.
    %
    % Input:
    %   thetas - nFx1 vector of angles relative to XY axis
    %   local_frames - 2*nFx3 matrix of local frames, first e1, then e2.
    %   N - nrosy field degree
    %   fids_list - if given, then the given angles are in these faces.
    %
    % Output:
    %   langle - angles relative to local frames
    
    n_angles = length(fids_list);
    theta = theta(:);
    %local_vectors = [cos(theta), sin(theta), zeros(n_angles, 1)];
    
    nF = size(local_frames, 1) / 2;
    if nargin < 4
        fids_list = 1:nF;
    end
    %ffield = zeros(N*n_angles, 3);
    langle = zeros(N*n_angles, 1);
    %for fid = fids_list
    for k = 1:n_angles
        fid = fids_list(k);
        % Assume frame is orthogonal
        frame = [local_frames(fid, :); ...
                 local_frames(fid+nF, :)];
        for i = 0:(N-1)
            t = theta(k) + i*2*pi/N;
            vec = frame * [cos(t), sin(t), 0]';
            %vec = frame * [cos(t); sin(t); 0];
            %ffield(k+i*n_angles, :) = vec;
            langle(k+i*n_angles) = atan2(vec(2), vec(1));
        end
    end

end

