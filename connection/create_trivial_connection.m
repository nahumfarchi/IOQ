function [x_star, E] = create_trivial_connection(m, ...
    vert_sing, ...
    gen_sing, ...
    fids_list, ...
    thetas_list, ...
    boundary_constraints, ...
    verbose)
    % function [x_star, E] = create_trivial_connection(m, S, fids_list, thetas_list, verbose)
    %
    % Computes the trivial connection as the solution to
    %   [d0'; H'] * x = -b
    % where 
    %   - d0 is the discrete exterior derivative on 0-forms (each column is a
    %   cycle around a vertex with +-1 according to orientation), 
    %   - H is a matrix composed of the holonomy cycles in the columns.
    %   - b_j = K_j - 2*pi/N * I_j
    %       K_j is the Gaussian curvature of vertx j
    %
    % Input:
    %   m - mesh object
    %   I - (V+2g) x 1 singularity indices. Must satisfy that the sum{I_j}=\xi
    %       where the sum is over vertices and boundary loops.
    %   I - (nSing) x 2 singularity indices. Each row is a pair
    %       (vid/holonomy_id, singularity_index). Must satisfy that sum{I_j}=\xi,
    %       where the sum is over vertices and boundary loops.
    %   N - nrosy field degree. Equal to 1 by default.
    %
    % Output: 
    %   x_star - E x 1 connection vector.
    %   E - energy = |x_star|_2^2, where x is the connection

    if nargin < 5
        verbose = false;
    end
    
    if verbose
        if abs(m.genus) > 1e-10
            warning(['Currently works with genus 0 only, but genus is ', num2str(m.genus)])
        end

        xi = 2 - 2*m.genus;
        %singularity_sum = sum(I(I(:,1) <= m.nV, 2))
        singularity_sum = sum(vert_sing(:, 2));
        if abs(singularity_sum - xi) > 1e-10
            warning(['Singularity indices must sum to ', num2str(xi), ' but sum to ', num2str(singularity_sum), '!'])
        end
    end   
   
    if verbose, disp('Creating connection...'); end
    
    [d0, d1] = get_exterior_derivatives(m);
       
    % Vertices
    b = get_gaussian_curvature(m);
    for j = 1:size(vert_sing,1)
        id = vert_sing(j, 1);
        b(id) = (b(id) - 2 * pi * vert_sing(j, 2));
    end
    
    % Get rid of vertices that are on the boundary
    EV = m.EVAdj;
    boundary_vert = zeros(m.nV, 1);
    be = m.boundaryEdge;
    for eid = 1:m.nE
        if be(eid)
            v1 = EV(eid, 1);
            v2 = EV(eid, 2);
            boundary_vert(v1) = true;
            boundary_vert(v2) = true;
        end
    end
    d0 = d0(:, ~boundary_vert);
    b = b(~boundary_vert);
    
    % Generators
    H = [];
    if abs(m.genus) > 1e-10
        H = m.H;
        %bg = generator_angle_defects(m);
        %bg = [1.9157; 2.2176];
        %bg = [0; 0];
        bg = -flatten_generators(m);
        %bg = [1.9877e-05; 1.5688e-05];
        %bg = [2.547155202214618; 2.01516679132591; 0.2044185474920219; -1.600589058173675];
        
        for j = 1:size(gen_sing,1)
            id = gen_sing(j, 1);
            bg(id) = (bg(id) - 2 * pi * gen_sing(j, 2));
        end
        b = [b; bg];
    end
    
    % System
    A = [d0'; -H'];
    
    % Constraints
    if ~isempty(boundary_constraints)
        A = [A; boundary_constraints.A];
        b = [b; boundary_constraints.b];
    end
    
    Z = d1;
    if length(fids_list) > 1
        [C, bc] = create_constraints_mat(m, fids_list, thetas_list);
        A = [A; C];
        b = [b; -bc];
        %Z = null(full(A))'; % There must be a more efficient way to do this...
        %fun = @(x) norm(x,2)^2;
        %x_tilda = A \ -b;
        %options = optimoptions('fmincon', 'MaxFunctionEvaluations', 7000);
        %x_star = fmincon(fun, x_tilda, [], [], A, -b, [], [], [], options);
    end
    if length(fids_list) > 1 || abs(m.genus) > 1e-10
        Z = nulls(A)';
    end
    
    % Project out null space component
%     x_tilda = A \ -b;
%     %%% Project out the null space component by subtracting 
%     %%%   Z' * inv(Z*Z') * Z * x_tilda
%     tmp = (Z * Z') \ (Z * x_tilda);
%     x_star = x_tilda - Z' * tmp;
%     check_norm('A*x_star', '-b');
    
    % QR decomposition
    % A' E = Q R
    % A x = -b
    % E R' Q x = -b
    % x = - (E R' Q)^{-1} b
    %[Q, R, E] = qr(A');
    %x_star = Q * (R' \ (E' * -b));
    %check_norm('A*x_star', '-b');
    
%     best_E = inf;
%     ng = size(H, 2);
%     for i = npermutek([1, -1], ng)'
%         b_tmp = b;
%         b(end:-1:end-ng+1) = b(end:-1:end-ng+1) .* i;
%         [Q, R, E] = qr(A');
%         x_star = Q * (R' \ (E' * -b));
%         check_norm('A*x_star', '-b');
%         norm(x_star)^2
%         if norm(x_star)^2 < best_E
%             best_E = norm(x_star)^2'
%             x_best = x_star;
%         end
%         b = b_tmp;
%     end
%     x_star = x_best;

    % least squares
    n = size(A, 2);
    C = speye(n, n);
    d = zeros(n, 1);
    global bk;
    -b'
    %b(1:end-4) = -b(1:end-4);
    [A, K, d0, d1, H] = tcods_gsystem(m.V, m.F);
    ng2 = size(H, 2);
    b = -K;
    b(vert_sing(:, 1)) = b(vert_sing(:, 1)) + 2 * pi * vert_sing(:, 2);
    b(m.nV + gen_sing(:, 1)) = b(m.nV + gen_sing(:, 1)) + 2 * pi * gen_sing(:, 2);
    x_star = lsqlin(C, d, [], [], A, b, -inf(n, 1), inf(n, 1));
    check_norm('A*x_star', '-b', 'Log', -1);
    
    
    if verbose
        %check_norm('A*x_tilda', '-b');
        %check_norm('A*x_star', '-b');
        %disp(['|| x_tilda ||^2 = ', num2str(norm(x_tilda)^2)])
        %disp(['|| x_star ||^2 = ', num2str(norm(x_star)^2)])
        %disp(['|| A * x_tilda - b ||^2 = ', num2str(norm(A * x_tilda - (-b))^2)])
        %disp(['|| A * x_star - b ||^2 = ', num2str(norm(A * x_star - (-b))^2)])
        %disp('Done.')
        
        
    end
    
    E = norm(x_star)^2;
end

