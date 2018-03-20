function [unwrapped, n, P, Q] = unwrap_phases_grid(V, E, wrapped, N, ub, lb)
    % function [unwrapped, n, edge_periods] = unwrap_phases(m, wrapped, frame_diffs, N, ub, lb)
    %
    % Phase unwrapping method based on the edgelist formulation [Piyush & Zebker].
    %
    % 
    % | 1 | 2 | ... | 2*nE      | 2*nE+1      | ... | 2*nE+nV      |
    % |------------------------------------------------------------|
    % |        P_ij, Q_ij       |                n_i               |
    % |------------------------------------------------------------|
    %
    % Q_ij - even indices
    % P_ij - odd indices
    %
    % Solve:
    %   min sum( cost_ij*(P_ij+Q_ij) )
    %   s.t.
    %       n_j-n_i+P_ij-Q_ij = round( (w_i-w_j)*N/(2pi) )
    %       n_i - integer
    %       P_ij,Q_ij >= 0 and real
    %   (for all edges (i,j))
    %
    % Input:
    %   V - grid vertices
    %   E - grid edges
    %   wrapped - wrapped phases
    %   N - directional field degree
    %   ub, lb - upper and lower bound constraints on flow
    %
    % Output:
    %   unwrapped - unwrapped phases
    %   n - n_i variables
    %   P, Q - integer flow
    
    nE = 0;
    for i = 1:size(E, 1)
        if(E(i, 1) > 0 && E(i,2) > 0)
            nE = nE + 1;
        end
    end
    
    %nE = size(E, 1);
    %n_rows = size(E, 1);
    n_rows = nE;
    n_cols = 2*n_rows + size(V, 1);
    
    cost = ones(n_cols, 1);
    lb = zeros(n_cols, 1);
    ub = 50*ones(n_cols, 1);
    b = zeros(n_rows, 1);
    
    I = [];
    J = [];
    S = [];
    offset = 2*n_rows;
    row = 1;
    
    for eid = 1:size(E, 1) %n_rows
        v1 = E(eid, 1);
        v2 = E(eid, 2);
        if v1 < 1 || v2 < 1
            continue
        end
        
        % n_i
        I(end+1) = row;
        J(end+1) = offset + v1;
        S(end+1) = -1;
        
        % n_j
        I(end+1) = row;
        J(end+1) = offset + v2;
        S(end+1) = 1;
        
        % P_ij
        col = 2*(row-1) + 1;
        I(end+1) = row;
        J(end+1) = col;
        S(end+1) = 1;
        
        % Q_ij
        col = 2*(row-1) + 2;
        I(end+1) = row;
        J(end+1) = col;
        S(end+1) = -1;
        
        % RHS
        w1 = wrapped(v1);
        w2 = wrapped(v2);
        b(row) = round((w1-w2)*N/(2*pi));
        
        row = row + 1;
    end
    
    A = sparse(I, J, S, n_rows, n_cols);
    
    [x, fval, exitflag] = linprog(cost, [], [], A, b, lb, ub);
    
    
    unwrapped = zeros(size(wrapped));
    for i = 1:length(wrapped)
        unwrapped(i) = wrapped(i)+2*pi*x(offset+i) / N;
    end
    
    if nargout > 1
        n = x(offset+1:end);
    end
    
    if nargout > 2
        P = x(2:2:2*nE);
    end
    
    if nargout > 3
        Q = x(1:2:2*nE);
    end

end

