function [theta] = representative_to_local(R, local_frames)
    %function [theta] = representative_to_local(R, local_frames)
    % Convert the given representative vectors (in 3D)
    % to angles in relation to the given local frames.
    %
    % Input:
    %   R - nF x 3 matrix of representative vectors
    %   local_frames - 2*nF x 3 matrix of local frames
    %
    % Output:
    %   theta - nF x 1 vector of angles with respect to the local frames
    
    nF = size(R, 1);
    theta = zeros(nF, 1);
    for i = 1:nF
        frame = [local_frames(i, :); ...
                 local_frames(i+nF, :)];
        v = frame * R(i,:)';
        theta(i) = atan2(v(2), v(1));
    end

end

