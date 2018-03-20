function [x, uk] = closest_lpoint_vor(B, y, u0)
    % function [x] = closest_lpoint_vor(B, y)
    %
    % Find the closest point in a lattice of Voronoi's first kind. Such a
    % lattice  has a set of vectors B = [b1,...,b_{n+1}] in R^n that satisfy
    %   1. b1,...,b_n form a basis for the lattice
    %   2. sum_{i=1}^{n+1} b_i = 0
    %   3. <b_i, b_j> <= 0 for i,j=1...n+1, i \neq j
    %
    % Implements https://arxiv.org/abs/1405.7014. Works in O(n^4).
    %
    % Input:
    %   B - nx(n+1) matrix, where the first n columns are the basis vectors
    %       of the lattice.
    %   y - n-dimensional point for which we want to find the closest 
    %       lattice point.
    
    EPS = 1e-10;
    
    n = size(B, 2) - 1;
    Q = B'*B;
    z = B(:,1:end-1)\y;
    % If z_i is zero then for numeric reasons it can either be +eps or
    % -eps, which is bad when rounding down.
    if nargin < 3
        u0 = floor(z+EPS); 
    end
    uk = u0;
    
    if any(any(Q-diag(diag(Q)) > 0))
        warning('Off diagonal elements of Q should be <= 0: %g', max(Q-diag(diag(Q))))
    end
    if any(sum(B, 2) > EPS)
        warning('Columns of B should sum to zero: %g', max(sum(B,2)));
    end
    
    for k = 1:n
        % Minimize ||B(z-uk-t)||_2 over t \in {0,1}^{n+1}
        t = solve_for_t(Q, z, uk);
        uk = uk + t;
    end
    
    x = B*uk;
    
function t = solve_for_t(Q, z, uk)
    % Solve the network flow problem
    %   w_ij      = -q_ij for i,j=2,...,n+2
    %   w_{i,n+3} = s(i-1) for i=2,...,n+2 and if s(i-1)>=0
    %               0 otherwise
    %   w_{1,i}   = 0 for i=2,...,n+2 and if s(i-1)>=0
    %             = -s(i-1) otherwise
    %
    % The minimum cut C gives t:
    %   t_i = 1 if i \in C
    %       = 0 otherwise
    
    n = size(Q, 2) - 1;
    W = zeros(n+3, n+3);
    
    p = z - uk;
    s = -2*Q*p;
    
    W(2:end-1, 2:end-1) = -Q;
    diag_inds = 1 : n+3+1 : (n+3)^2;
    W(diag_inds(2:end-1)) = diag(Q);
    
    W(1, 1) = 0;
    W(end, end) = 0;
    W(1, end) = 0;
    W(end, 1) = 0;
    for i = 2:n+2
        if s(i-1) >= 0
            W(i, n+3) = s(i-1);
            W(1, i) = 0;
        else
            W(i, n+3) = 0;
            W(1, i) = -s(i-1);
        end
    end
    W(:, 1) = W(1, :);
    W(end, :) = W(:, end);
    
    G = graph(abs(W));
    [~,~,cs,~] = maxflow(G,1,n+3);
    t = zeros(n+3, 1);
    t(cs) = 1;
    t = t(2:end-1);
    
function res = Qf(Q, s, t, n)
    res = 0;
    for i_ = 1:n+1
        res = res + s(i_)*t(i_);
    end
    for i_ = 1:n+1
        for j_ = 1:n+1
            res = res+Q(i_,j_)*t(i_)*t(j_);
        end
    end
