function [] = exp_2017_08_31_1145_iterative()

    global LOG;
    global EPS;
    global VERBOSE;

    ME = [];
    try
        ABOUT = 'VF using iterative heuristic\r\n\r\n';
        OUT_FOLDER_NAME = 'iterative';

        %MESHES = {'sphere_s0.off', ...
        %    'bumpy.off', ...
        %    'cow.off', ...
        %    'bunny.off', ...
        %    'bunny2.off', ...
        %    'phands.off', ...
        %    'torus_fat_r2.off', ...
        %    'round_cuber.off', ...
        %    'rounded_cube_keenan.off'};

        MESHES = {'rounded_cube_keenan.off'};

        %VIEW_ANGLE = [];

        VERBOSE = true;
        EPS = 1e-9;

        PLOT = true;
        SAVE = true;
        N_ITER = 100;

        if PLOT && SAVE
            OUT_FOLDER = create_time_stamped_folder(...
                fullfile('..', 'results'), ...
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

        DEGREE = 4;
        F0 = [1]; % Constrained face
        THETA0 = 0; % Constraint angle

        for fname = MESHES
            log_and_print(LOG, 'Loading %s...\r\n', fname{:});
            p = find_data_folder();
            fp = fullfile(p, fname{:});

            mesh = Mesh();
            mesh.loadTM(fp);
            V0 = mesh.V(mesh.F(F0, 2), :) - mesh.V(mesh.F(F0, 1), :);   

            log_and_print(LOG, 'degree : %d\r\n', DEGREE);

            res_MIQ = run_MIQ(fp, F0, V0, DEGREE);
            %res_GO = run_GO(fp, F0, V0, DEGREE);
            %res_TC_GO = TCODS(mesh, res_GO.S, F0, THETA0, DEGREE, VERBOSE);
            res_iterative = run_lattice_iterative(mesh, ...
                F0, ...
                THETA0, ...
                DEGREE, ...
                N_ITER);
            
            %res_TC_GO.title = {'TC (w/ GO singularities)', ...
            %    sprintf('$E_{MIQ} = %g$', res_TC_GO.E)};

            %results = {res_MIQ, res_iterative, res_GO, res_TC_GO};
            %results = {res_iterative};
            results = {res_MIQ, res_iterative};
            
            MeshVis.wfigs(fname{:}, ...
                        mesh, ...
                        'OutFolder', OUT_FOLDER, ...
                        'Nrosy', results, ...
                        'ConstrainedFaces', F0, ...
                        'ConstraintVectors', V0, ...
                        'Save', SAVE);
        end

    % A hack since matlab does not have a finally clause
    catch ME
    end
    % Close open resources
    if LOG > 0
        fclose(LOG);
    end
    if ~isempty(ME)
        rethrow(ME);
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

function res = run_MIQ(fp, f0, v0, degree)
    global LOG;
    log_and_print(LOG, '\r\nSolving with libigl MIQ...\r\n');
    res = nrosy_wrapper(fp, f0, v0, degree);
    log_and_print(LOG, '%s\r\n', res.result);
    log_and_print(LOG, 'Status: %d\r\n', res.status);
    log_and_print(LOG, 'E: %g\r\n', res.E);
    
    elapsed = regexp(res.result, 'Total\s+(\d+(?:\.\d+)?)', 'match');
    elapsed = regexp(elapsed, '(\d+(?:\.\d+)?)', 'match');
    
    res.title = {'MIQ ', ...
        sprintf('$E_{MIQ} = %g$', res.E), ...
        sprintf('Elapsed time: %s', elapsed{:}{:})};
    
function res = run_GO(fp, f0, v0, degree)
    tic
    
    [GO_energy, MIQ_energy, x, ffield, S] = calc_GO_and_MIQ_energies(fp, degree, f0, v0);
    
    res.E_GO = GO_energy;
    res.E_MIQ = MIQ_energy;
    
    res.ffield = ffield;
    res.degree = degree;
    res.S = S;
    
    elapsed = toc;
    
    res.title = {'GO ', ...
        sprintf('$E_{MIQ} = %g$', full(res.E_MIQ)), ...
        sprintf('$E_{GO} = %g$', full(res.E_GO)), ...
        sprintf('Elapsed time: %s', elapsed)};
    
function res = run_lattice_iterative(mesh, f0, theta0, degree, n_iter)
    % see notes/2017-08-31 Iterative impr[3852].pdf
    tic
    global EPS;
    global VERBOSE;
    
    nV = mesh.nV;
    [d0, ~] = get_exterior_derivatives(mesh);
    Ad = get_gaussian_curvature(mesh);
    xi = 2 - 2*mesh.genus;
    
    L = d0'*d0;
    
    tic
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
    Lp = V*Dinv*V';
    toc
    
    % faster way to get the pseudo inverse   
    tic; Lp = inv(L + 1/nV) - 1/nV; toc
    tic; Lgpu = gpuArray(L+1/nV);
    Lp = inv(Lgpu) - 1/nV; toc

    k0 = (2/pi)*Ad;
    c = 4*xi;
    
    k = zeros(nV, 1);
    inds = randperm(nV, c);
    k(inds) = 1;
    
    % Resistence distance matrix
    Ldii = repmat(diag(Lp), 1, nV);
    R = Ldii + Ldii' - 2*Lp;
    b = 2*Lp*(k-k0);
    
    for iter = 1:n_iter
        
        %tic
        %[X, Y] = meshgrid(b, b);
        %D1 = (X' - Y');
        %toc
        %tic
        %D = b - b';
        %toc
        %assert(check_norm('D1', 'D'));
        D = bsxfun(@minus, b, b');
    
        Z = D + R;
        [i, j] = find(Z == min(Z(:)), 1);
    
        k_tilde = k;
        k_tilde(i) = k_tilde(i) + 1;
        k_tilde(j) = k_tilde(j) - 1;
    
        %E = (k-k0)'*Lp*(k-k0);
        %E_tilde = (k_tilde-k0)'*Lp*(k_tilde-k0);
        
        %if abs(E - E_tilde) < EPS
        %    break
        %end
        %assert(check_norm('abs(E-E_tilde)', 'abs(b(i)-b(j)+R(i,j))'));
        if b(i) - b(j) + R(i, j) > 0
            %assert(abs(E - E_tilde) < EPS)
            break;
        end
        
        k = k_tilde;
        b = b + 2*Lp(:,i) - 2*Lp(:,j);
    end
    
    inds = find(abs(k) > EPS);
    S = [inds, k(inds) / degree];
    res = TCODS(mesh, S, f0, theta0, degree, VERBOSE);
    
    fprintf('E = %g\n', res.E);
    
    elapsed = toc;
    res.title = {'Iterative', ...
        sprintf('$E_{MIQ} = %g$', res.E), ...
        sprintf('Elapsed time: %g', elapsed)};
