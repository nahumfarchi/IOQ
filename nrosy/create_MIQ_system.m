function [A, b, thetas, periods, theta_tags, period_tags] ...
    = create_MIQ_system(...
        mesh, ...
        N, ...
        constrained_faces, ...
        constraint_vectors, ...
        local_frames, ...
        frame_diffs, ...
        periods, ...
        period_is_fixed)
%function [A, b, thetas, periods, theta_tags, period_tags] = create_MIQ_system(mesh, N, constrained_faces, constraint_vectors, local_frames, frame_diffs, periods, period_is_fixed)
% Create MIQ system. If 'periods' and 'period_is_fixed'
% are not given then the search space is reduced as described in the MIQ
% paper. 
%
%   E_MIQ = ||Ax - b||^2
%
% Input:
%   mesh
%   N - NRosy degree
%   constrained_faces - nC x 1 vector of constraint ids
%   constraint_vectors - nC x 3 matrix of constraint vectors
%   local_frames - (2*nF) x 3 matrix of local frames
%   frame_diffs - nE x 1 vector of frame differences between adjacent faces
%   periods (optional) - nE x 1 vector of period jumps
%   period_is_fixed (optional) - nE x 1 boolean vector
%
% Output:
%   A - nE x (n_free_theta +
%   n_free_periods) sparse matrix
%   b  - nE x 1 vector
%   theta_tags - maps fids to position in the solution vector x, i.e.,
%       theta(fid) = x(theta_tags(fid))
%   period_tags - maps eids to position in solution vector x, i.e.,
%       periods(fid) = x(period_tags(eid))

    if nargin < 5
        [local_frames, frame_diffs] = create_local_frames(mesh);
    end
    if nargin < 7
        reduce_space = true;
    else
        reduce_space = false;
    end

    nC = length(constrained_faces);
    nF = mesh.nF;
    nE = mesh.nE;
    EV = mesh.EVAdj;
    EF = mesh.EFAdj;
    assert(size(EV,1) == nE);

    % Convert constraints to local angles
    thetas = zeros(nF, 1);
    face_is_constrained = zeros(nF, 1);
    for i = 1:nC
        fid = constrained_faces(i);
        v = constraint_vectors(i, :);
        frame = [local_frames(fid, :);
                 local_frames(fid+nF, :)];
        v_local = frame * v';
        thetas(fid) = atan2(v_local(2), v_local(1));
        face_is_constrained(fid) = true;
    end

    % Count and tag variables
    n_thetas = nF - nC;
    count = 1;
    theta_tags = -ones(nF, 1);
    for fid = 1:nF
        if ~face_is_constrained(fid)
            theta_tags(fid) = count;
            count = count + 1;
        end
    end
    assert(n_thetas == count - 1);

    % Reduce search space
    period_tags = -ones(nE, 1);
    if reduce_space
        [period_is_fixed, periods] = reduce_search_space(...
                        mesh, ...
                        face_is_constrained, ...
                        thetas, ...
                        frame_diffs);
    end

    n_periods = length(period_is_fixed) - sum(period_is_fixed);
    for eid = 1:nE
        if ~period_is_fixed(eid)
            period_tags(eid) = count;
            count = count + 1;
        end
    end

    disp(['n_period : ', num2str(n_periods)])
    disp(['n_theta : ', num2str(n_thetas)])
    
    n_variables = n_thetas + n_periods;
    assert(n_variables == count - 1);
    b = zeros(nE, 1);
    I = [];
    J = [];
    V = [];
    row = 0;
    
    for eid = 1:nE
        i = EF(eid, 1);
        j = EF(eid, 2);
        if mesh.isBoundaryEdge(eid) || face_is_constrained(i) && face_is_constrained(j)
            continue
        end
        
        row = row + 1;
        isFixed_i = face_is_constrained(i);
        isFixed_j = face_is_constrained(j);
        isFixed_p = period_is_fixed(eid);
        ti = thetas(i);
        tj = thetas(j);
        kij = frame_diffs(eid);
        pij = periods(eid);
        
        if ~isFixed_i
            I(end+1) = row;
            J(end+1) = theta_tags(i);
            V(end+1) = 1;
        else
            b(row) = b(row) - ti;
        end
        
        if ~isFixed_j
            I(end+1) = row;
            J(end+1) = theta_tags(j);
            V(end+1) = -1;
        else
            b(row) = b(row) + tj;
        end
        
        if ~isFixed_p
            I(end+1) = row;
            J(end+1) = period_tags(eid);
            V(end+1) = 2*pi / N;
        else
            b(row) = b(row) - 2*pi*pij / N;
        end
        
        b(row) = b(row) - kij;
        
    end
    
    %assert(row == n_variable);
    A = sparse(I, J, V, nE, n_variables);

end

