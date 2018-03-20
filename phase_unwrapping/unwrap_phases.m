function [unwrapped, n, edge_periods] = unwrap_phases(m, wrapped, frame_diffs, N, ub, lb)
    % function [unwrapped, n, edge_periods] = unwrap_phases(m, wrapped, frame_diffs, N, ub, lb)
    %
    % Phase unwrapping method based on the edgelist formulation.
    %
    % 
    % | 1 | 2 | ... | 2*nE_dual | 2*nE_dual+1 | ... | 2*nE_dual+nF |
    % |------------------------------------------------------------|
    % |        P_ij, Q_ij       |                n_i               |
    % |------------------------------------------------------------|
    %
    % Q_ij - at even indices
    % P_ij - at odd indices
    %
    %
    % Input:
    %   m - mesh object
    %   wrapped - wrapped phases
    %   frame_diffs - angle differences between adjacent frames
    %   N - directional field degree
    %   ub, lb - upper and lower bound constraints on flow
    %
    % Output:
    %   unwrapped - unwrapped phases
    %   n - n_i variables
    %   edge_periods - for each edge, the difference between n_i and n_j
    
    if nargin > 4
        error('ub, lb are a mess')
    end
    
    %m = self.mesh;
    nE_dual = m.niE;
    nV_dual = m.nF;
    n_rows = nE_dual;
    n_cols = 2*nE_dual + nV_dual;
    offset = 2*nE_dual;
    %wrapped = self.thetas;

    EF = m.EFAdj;

    I = [];
    J = [];
    V = [];
    cost = zeros(n_cols, 1);
    lb = zeros(n_cols, 1);
    ub = zeros(n_cols, 1);
    beq = zeros(n_rows, 1);
    row = 1;
    for eid = 1:m.nE
        if m.isBoundaryEdge(eid)
            continue;
        end

        %fid1 = min(EF(eid, 1), EF(eid, 2));
        %fid2 = max(EF(eid, 1), EF(eid, 2));
        fid1 = EF(eid, 1);
        fid2 = EF(eid, 2);

        % n_i
        I(end+1) = row;
        J(end+1) = offset + fid1;
        V(end+1) = -1;
        % n_j
        I(end+1) = row;
        J(end+1) = offset + fid2;
        V(end+1) = 1;
        % P_ij
        col = 2*(row-1) + 1;
        I(end+1) = row;
        J(end+1) = col;
        V(end+1) = 1;
        cost(col) = 1;
        lb(col) = 0;
        ub(col) = 50;
        % Q_ij
        col = 2*(row-1) + 2;
        I(end+1) = row;
        J(end+1) = col;
        V(end+1) = -1;
        cost(col) = 1;
        lb(col) = 0;
        ub(col) = 50;      
        % RHS
        t1 = wrapped(fid1);
        t2 = wrapped(fid2);
        rij = frame_diffs(eid);
        %period_ij = self.periods(eid);
        beq(row) = round((t1-t2+rij) * N / (2*pi));
        %beq(row) = round((t1-t2)/(2*pi));

        row = row + 1;
    end

    cost(offset+1:end) = 0;
    lb(offset+1:end) = -50;
    ub(offset+1:end) = 50;

    Aeq = (sparse(I, J, V, n_rows, n_cols));

    [x,fval,exitflag] = linprog(cost, [], [], Aeq, beq, lb, ub);
    x = round(x, 0);
    
    unwrapped = zeros(size(wrapped));
    for fid = 1:m.nF
        unwrapped(fid) = wrapped(fid) + 2*pi*x(offset+fid) / N;
    end

    n = x(offset+1:end);
    
    edge_periods = zeros(m.nE, 1);
    for eid = 1:m.nE
        fid1 = EF(eid, 1);
        fid2 = EF(eid, 2);
        %self.periods(eid) = self.periods(eid) + (x(offset+fid1) - x(offset+fid2));
        edge_periods(eid) = (x(offset+fid1) - x(offset+fid2));
    end      

end

