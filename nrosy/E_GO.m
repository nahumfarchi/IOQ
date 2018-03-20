function [E] = E_GO(nrosy, local_frames)
    % function [E] = E_GO(nrosy)
    % Calculate GO energy.
    %
    % Input:
    %   struct nrosy with the fields:
    %       theta - nFx1 vector of angles
    %       ffield - (nF*degree)x3 vector field on the faces
    nF = length(nrosy.theta);
    u = langles_to_complex(nrosy.theta);
    for fid = 1:nF
        v = nrosy.ffield(fid, :);
        c = dot(v, local_frames(fid,:)) + 1i*dot(v, local_frames(fid+mesh.nF,:));
        u(fid) = c;
    end
    E = norm(A*u-b)^2;

end

