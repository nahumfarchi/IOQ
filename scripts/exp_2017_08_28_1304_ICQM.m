function exp_2017_08_28_1304_ICQM()
    % Example script for running ICQM

    % Construct integer least squares problem
    % minimize    ||Ax - b||^2
    % subject to  x in Z^n

    %n = 100;
    %A = randn(2*n, n) / 10; b = randn(2*n, 1);

    ABOUT = 'Lattice using ICQM method\r\n\r\n';
    OUT_FOLDER_NAME = 'lattice_ICQM';

    %MESHES = {'sphere_s0.off', ...
    %    'bumpy.off', ...
    %    'round_cuber.off', ...
    %    'rounded_cube_keenan.off', ...
    %    'cow.off', ...
    %    'bunny.off', ...
    %    'bunny2.off', ...
    %    'phands.off', ...
    %    'torus_fat_r2.off'};

    %MESHES = {'bunny2.off', 'bunny.off', 'rounded_cube_keenan.off', 'round_cuber.off'};
    %MESHES = {'rounded_cube_keenan.off'};
    %MESHES = {'round_cuber.off'};
    %MESHES = {'sphere_s0.off'};
    %MESHES = {'sphere_s0.off'};
    MESHES = {'bumpy.off'};

    %VIEW_ANGLE = [];

    VERBOSE = true;
    EPS = 1e-9;

    theta0 = 0;

    PLOT = true;
    SAVE = true;

    if PLOT && SAVE
        OUT_FOLDER = create_time_stamped_folder(fullfile('..', 'results'), ...
            OUT_FOLDER_NAME, ...
            true, ...
            false);
        LOG = fopen(fullfile(OUT_FOLDER, 'log.txt'), 'w');
    else
        LOG = -1;
    end

    log_and_print(LOG, '%s\r\n', mfilename('fullpath'));
    log_and_print(LOG, ABOUT);
    log_and_print(LOG, 'logid : %d\r\n', LOG);

    for fname = MESHES
        log_and_print(LOG, 'Loading %s...\r\n', fname{:});
        p = find_data_folder();
        fp = fullfile(p, fname{:});

        mesh = Mesh();
        mesh.loadTM(fp);

        f0 = [1]; % Starting face
        v0 = mesh.V(mesh.F(f0, 2), :) - mesh.V(mesh.F(f0, 1), :);
        degree = 4;

        log_and_print(LOG, '\r\nSolving with libigl MIQ...\r\n');
        res_MIQ = nrosy_wrapper(fp, f0, v0, degree);
        log_and_print(LOG, '%s\r\n', res_MIQ.result);
        log_and_print(LOG, 'Status: %d\r\n', res_MIQ.status);
        log_and_print(LOG, 'E: %g\r\n', res_MIQ.E);

    %     %
    %     %
    %     %
    %     [A, b, thetas, periods, theta_tags, period_tags] = ...
    %     create_MIQ_system(mesh, 4, [1], [1,0,0]);
    % 
    %     % [A11 A12] [theta] | [b1]
    %     % [A21 A22] [p    ] | [b2]
    %     theta_inds = theta_tags(theta_tags>0);
    %     period_inds = period_tags(period_tags>0);
    %     A11 = A(theta_inds, theta_inds);
    %     A12 = A(theta_inds, period_inds);
    %     A21 = A(period_inds, theta_inds);
    %     A22 = A(period_inds, period_inds);
    %     b1 = b(theta_inds);
    %     b2 = b(period_inds);
    %     
    %     S = A22 - A21 * (A11 \ A12);
    %     b_tilde = b2 - A21 * (A11 \ b1);
    %     
    %     % ||Ax - b||^2 = x^T P x + 2q^T x + r
    %     P = S'*S;
    %     q = -S*b_tilde;
    %     r = b_tilde'*b_tilde;
    %     
    %     tic;
    %     [lb, ub, xhat] = ICQM(P, q, r);
    %     elapsed = toc;
    %     
    %     fprintf('Elapsed time: %.5f seconds\n', elapsed);
    %     fprintf('Optimal value is between %.5f and %.5f\n', lb, ub);
    %     fprintf('2-norm of the suboptimal solution: %.5f\n', norm(xhat));

        %
        %
        %
        nV = mesh.nV;
        [d0, ~] = get_exterior_derivatives(mesh);
        Ad = get_gaussian_curvature(mesh);
        L = d0'*d0;
        [V,D] = eig(full(L));
        [D, sort_eigen] = sort(diag(D));
        V = V(:, sort_eigen);
        D = diag(D);

        D = D(1:nV, 1:nV);
        V = V(:, 1:nV);

        Dinv = diag(D);
        Dinv(Dinv < EPS) = inf;
        Dinv = 1 ./ Dinv;
        Dinv = diag(Dinv);
        Ldagger = V*Dinv*V';

        % |M*x-y|
        % |M*(x-z)|
        M = V*sqrt(Dinv)*V';
        %M = d0*Wpinv;
        %b = (2/pi)*Ad;
        y = (2/pi)*M*Ad;
        %m = 4*xi;

        xi = 2 - 2*mesh.genus;
        c = 4*xi;
        n = size(M, 2);
        [M, y] = add_constraint(M, y, c);

        P = M'*M;
        q = -M'*y;
        r = y'*y;

        tic;
        [lb, ub, xhat] = ICQM(P, q, r);
        elapsed = toc;

        k = [xhat; c - ones(1, n-1)*xhat];
        S = find(abs(k)>EPS);
        S = [S, k(S)/degree];
        res_lattice_ICQM = TCODS(mesh, S, [1], 0, degree, VERBOSE);

        log_and_print(LOG, 'Elapsed time: %.5f seconds\n', elapsed);
        log_and_print(LOG, 'Optimal value is between %.5f and %.5f\n', lb, ub);
        log_and_print(LOG, '2-norm of the suboptimal solution: %.5f\n', norm(xhat));

        %
        %
        %
        nV = mesh.nV;
        [d0, ~] = get_exterior_derivatives(mesh);
        Ad = get_gaussian_curvature(mesh);
        L = d0'*d0;
        [V,D] = eig(full(L));
        [D, sort_eigen] = sort(diag(D));
        V = V(:, sort_eigen);
        D = diag(D);

        D = D(1:nV, 1:nV);
        V = V(:, 1:nV);

        Dinv = diag(D);
        Dinv(Dinv < EPS) = inf;
        Dinv = 1 ./ Dinv;
        Dinv = diag(Dinv);
        Ldagger = V*Dinv*V';

        % |M*x-y|
        % |M*(x-z)|
        %M = V*sqrt(Dinv)*V';
        M = d0*Ldagger;
        %b = (2/pi)*Ad;
        y = (2/pi)*M*Ad;
        %m = 4*xi;

        xi = 2 - 2*mesh.genus;
        c = 4*xi;
        n = size(M, 2);
        [M, y] = add_constraint(M, y, c);

        % ||Mx - y||^2 = x^T P x + 2q^T x + r
        P = M'*M;
        q = -M'*y;
        r = y'*y;

        tic;
        [lb, ub, xhat] = ICQM(P, q, r);
        elapsed = toc;

        k = [xhat; c - ones(1, n-1)*xhat];
        S = find(abs(k)>EPS);
        S = [S, k(S)/degree];
        res_lattice_ICQM2 = TCODS(mesh, S, [1], 0, degree, VERBOSE);

        log_and_print(LOG, 'Elapsed time: %.5f seconds\n', elapsed);
        log_and_print(LOG, 'Optimal value is between %.5f and %.5f\n', lb, ub);
        log_and_print(LOG, '2-norm of the suboptimal solution: %.5f\n', norm(xhat));

        res_MIQ.title = {'MIQ (libigl)', sprintf('$E_{MIQ} = %g$', res_MIQ.E)};
        res_lattice_ICQM.title = {'$\|V\sqrt{D^{-1}}V^{T}x - (2/pi)V\sqrt{D^{-1}}V^{T}A_d\|^2$ solved with ICQM', ...
            sprintf('$E_{MIQ} = %g$', res_lattice_ICQM.E)};
        res_lattice_ICQM2.title = {'$\|d_0(d_0^Td_0)^{\dagger}x - (2/pi)d_0(d_0^Td_0)^{\dagger}A_d\|^2$ solved with ICQM', ...
            sprintf('$E_{MIQ} = %g$', res_lattice_ICQM2.E)};

        MeshVis.wfigs(fname{:}, ...
                mesh, ...
                'OutFolder', OUT_FOLDER, ...
                'Nrosy', {res_MIQ, res_lattice_ICQM, res_lattice_ICQM2}, ...
                'ConstrainedFaces', f0, ...
                'ConstraintVectors', v0, ...
                'Save', SAVE);
    end
    
function [X, b] = add_constraint(M, y, c)
    % x_I^{*} = argmin |X x - b|
    %           x \in Z^{n-1}
    % x^{*} = [x_I^{*}; ones(1, n-1) x_I^{*}]
    n = size(M, 2);
    M_B = M(:, 1);
    M_I = M(:, 2:end);
    X = M_I - M_B*ones(1,n-1);
    b = y - M_B*c;
