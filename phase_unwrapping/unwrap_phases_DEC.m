function [h] = unwrap_phases_DEC(m, g, h1)
    %function [] = unwrap_phases_DEC(m, g)
    % Unwrap the phases on the given mesh. Assumes that m is flat.
    %
    % Input:
    %   m - the mesh
    %   g - wrapped phases
    %   h1 - start with h(v1) = h1. Default value is g(v1).
    % Output:
    %   h - unwrapped phases
    %
    % Main idea:
    % \alpha = d0*g + 2\pi\eta
    %
    % Minimize ||\beta||_1,
    % s.t. 
    %   d1\beta = -d1\eta
    %
    % Estimate dh by
    %   d0*h = 2\pi\beta+\alpha
    % and use it to reconstruct the unwrapped phases.

    if nargin < 3
        h1 = g(1); % Assume h(v1) = g(1)
    end

    [d0, d1] = get_exterior_derivatives(m);
    nE = m.nE;

    % Calculate alpha=d0*g+2pi\eta
    dg = d0*g;

    % find 2*pi multipels such that alpha \in [-pi,pi)
    tmp = (dg+sign(dg)*pi) / (2*pi);
    eta = zeros(size(tmp));
    eta(tmp>0) = -floor(tmp(tmp>0));
    eta(tmp<0) = -ceil (tmp(tmp<0));

    alpha = dg + 2*pi*eta;
    assert(~(any(alpha < -pi | alpha > pi)))

    % Solve the minimization problem
    cost = ones(2*nE, 1);
    lb = zeros(2*nE, 1);
    ub = 50*ones(2*nE, 1);
    A = [d1, -d1]; % \beta = xplus - xminus
    b = -d1*eta;
    [x, fval, exitflag] = linprog(cost, [], [], A, b, lb, ub);

    % Set beta to the non-zero elements of xplus and xminus
    xplus = x(1:nE);
    xminus = x(nE+1:end);
    beta = zeros(m.nE, 1);
    beta(abs(xplus) > 1e-10) = xplus(abs(xplus) > 1e-10);
    beta(abs(xminus) > 1e-10) = xminus(abs(xminus) > 1e-10);
    
    % To calculate the unwrapped phases,
    % run a bfs starting at v1, and for each edge e=(vi,vj) 
    % add dh(e) to vj
    dh = 2*pi*beta + alpha;

    VV_to_dh = sparse(m.nV, m.nV);
    EV = m.EVAdj;

    for e = 1:m.nE
        vi = EV(e, 1);
        vj = EV(e, 2);
        VV_to_dh(vi, vj) = dh(e);
        VV_to_dh(vj, vi) = -dh(e);
    end

    % Create stack
    import java.util.LinkedList
    q = LinkedList();
    q.addLast(1);
    visited = zeros(m.nV, 1);
    
    % Init unwrapped phases
    h = zeros(m.nV, 1);
    h(1) = h1;
    path_calculated = zeros(m.nV, 1);
    path_calculated(1) = true;

    % Sum up estimated grads by bfs
    VV = m.VVAdj;
    while ~q.isEmpty()
        vi = q.removeFirst();
        visited(vi) = true;

        for vj = find(VV(vi,:) ~= 0)
            if visited(vj)
                continue;
            end
            % If there was no path from v1 to vj yet, then update h(vj)
            if ~path_calculated(vj)
                h(vj) = h(vi) + VV_to_dh(vi, vj);
                path_calculated(vj) = true;
            end
            q.addLast(vj);
        end
    end
end

