function [] = exp_2017_08_31_1145_iterative()

    global LOG;
    global EPS;
    global VERBOSE;
    global PLOT;
    global SAVE;
    global DEGREE;
    global ASPECT_RATIO;
    global RESOLUTION;
    global OUT_FOLDER;
    VERBOSE = true;
    EPS = 1e-9;
    PLOT = true;
    SAVE = true;
    DEGREE = 4;
    ASPECT_RATIO = 1;
    RESOLUTION = 1024;
    PROFILE = false;
    
    OUT_FOLDER_NAME = 'iterative';
    
    if PROFILE
        profile clear;
        profile on;
    end
    
    ME = [];
    try
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
        
        log_and_print(LOG, 'degree : %d\r\n', DEGREE);
    
        %MESHES_ALL = {'sphere_s0.off', ...
            %    'bumpy.off', ...
            %    'cow.off', ...
            %    'bunny.off', ...
            %    'bunny2.off', ...
            %    'phands.off', ...
            %    'torus_fat_r2.off', ...
            %    'round_cuber.off', ...
            %    'rounded_cube_keenan.off'};

        %MESHES_CUBES = {'rounded_cube_keenan.off'};

        %fp = '../data/bunny_res7.off';
        %n_sbd = 3;
        %subdivide_and_save(fp, n_sbd);

        %timing_plot({'sphere_s0.off', 'sphere_s1.off', 'sphere_s2.off', 'sphere_s3.off'});
        %timing_plot({'sphere_s0.off', 'sphere_s1.off', 'sphere_s2.off', 'sphere_s3.off', 'sphere_s4.off'});
        %timing_plot({'sphere_s0.off', 'sphere_s1.off', 'sphere_s2.off', 'sphere_s3.off', 'sphere_s4.off', 'sphere_s5.off'});
        %timing_plot({'bunny_s0.off', 'bunny_s1.off'});
        %timing_plot({'bunny_simple_s0.off', 'bunny_simple_s1.off', 'bunny_simple_s2.off', 'bunny_simple_s3.off'});
        %energy_and_time_plots({'bunny_res1.off', 'bunny_res2.off', 'bunny_res4.off', 'bunny_res5.off', 'bunny_res6.off', 'bunny_res7.off'});
        %energy_and_time_plots({'bunny_res6.off', 'bunny_res7.off'});
        
        %profile_mesh('../data/bunny_res7_sd1.off');

        %create_and_save_plots(MESHES);
        
        %lap_pseudo_inverse_timings({'bunny_res1.off', 'bunny_res2.off', 'bunny_res3.off', 'bunny_res4.off', 'bunny_res5.off', 'bunny_res6.off', 'bunny_res7.off', 'bunny_res7_sd1.off'});
        %lap_pseudo_inverse_timings({'bunny_res7_sd1.off'});
        
        %timing_plot_MIQ_iter_eigs({'bunny_res1.off', 'bunny_res2.off', 'bunny_res3.off', 'bunny_res4.off', 'bunny_res5.off', 'bunny_res6.off', 'bunny_res7.off', 'bunny_res7_sd1.off'}, 'bunny');
        %timing_plot_MIQ_iter_eigs({'bunny_res1.off', 'bunny_res2.off', 'bunny_res3.off'}, 'bunny');
        %timing_plot_MIQ_iter_eigs({'bunny_res1.off'}, 'bunny');
        
        %gpu_memory_tests('../data/bunny_res7_sd2.off');
        %gpu_memory_tests('../data/bunny_57k_faces.off');
        
        test_CUDA('../data/bunny_res6.off');
    
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
    
    if PROFILE
        profile viewer;
    end
    
function profile_mesh(fp)
    global LOG;   
    global DEGREE;
    
    N_ITER = 1000;
    F0 = [1];
    THETA0 = [0];

    mesh = Mesh();
    mesh.loadTM(fp);   
    V0 = mesh.V(mesh.F(F0, 2), :) - mesh.V(mesh.F(F0, 1), :);   

    res_MIQ = run_MIQ(fp, F0, V0, DEGREE);
    %res_GO = run_GO(fp, F0, V0, DEGREE);
    %res_TC_GO = TCODS(mesh, res_GO.S, F0, THETA0, DEGREE, VERBOSE);

    opt.plot = true;
    opt.cot = true;
    opt.graph = false;

    opt.gpu = true;
    opt.chol = false;
    %profile clear; profile on
    res_iter_gpu_nochol = run_lattice_iterative(mesh, ...
        F0, ...
        THETA0, ...
        DEGREE, ...
        N_ITER, ...
        opt);
    %profile viewer
    
    log_and_print(LOG, 'Elapsed MIQ  : %g\r\n', res_MIQ.elapsed);
    log_and_print(LOG, 'Elapsed iter : %g\r\n', res_iter_gpu_nochol.elapsed_iter);
        
function create_and_save_plots(meshes)
    global LOG;
    global DEGREE;
    
    ABOUT = 'VF using iterative heuristic\r\n\r\n';
    %VIEW_ANGLE = [];
    N_ITER = 100;

    log_and_print(LOG, '%s\r\n', mfilename('fullpath'));
    log_and_print(LOG, ABOUT);
    log_and_print(LOG, 'logid : %d\r\n', LOG);

    F0 = [1]; % Constrained face
    THETA0 = 0; % Constraint angle

    for fname = meshes
        log_and_print(LOG, 'Loading %s...\r\n', fname{:});
        p = find_data_folder();
        fp = fullfile(p, fname{:});

        mesh = Mesh();
        mesh.loadTM(fp);
        V0 = mesh.V(mesh.F(F0, 2), :) - mesh.V(mesh.F(F0, 1), :);   

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

function timing_plot_MIQ_iter_eigs(mesh_names, out_filename)
    global LOG;   
    global DEGREE;
    global OUT_FOLDER;
    global SAVE;
    
    N_ITER = 1000;
    F0 = [1];
    THETA0 = [0];
    
    n_meshes = length(mesh_names);
    n_faces = zeros(n_meshes, 1);
    T_MIQ = zeros(n_meshes, 1);
    T_iter_gpu = zeros(n_meshes, 1);
    %T_eigs = zeros(n_meshes, 1);
    
    E_MIQ = zeros(n_meshes, 1);
    E_iter_gpu = zeros(n_meshes, 1);
    %E_eigs = zeros(n_meshes, 1);
    
    n_sing_MIQ = zeros(n_meshes, 1);
    n_sing_iter = zeros(n_meshes, 1);
    %n_sing_eigs = zeros(n_meshes, 1);
    
    results_MIQ = {};
    results_iter = {};
    %results_eigs = {};
    meshes = {};
    
    k = 1;
    fig1 = figure();
    fig2 = figure();
    fig3 = figure();
    for fname = mesh_names
        log_and_print(LOG, 'Loading %s...\r\n', fname{:});
        
        p = find_data_folder();
        fp = fullfile(p, fname{:});

        mesh = Mesh();
        mesh.loadTM(fp);   
        V0 = mesh.V(mesh.F(F0, 2), :) - mesh.V(mesh.F(F0, 1), :);   

        res_MIQ = run_MIQ(fp, F0, V0, DEGREE);
        res_MIQ.title = {res_MIQ.title{:}, ...
            sprintf('n sing: %d', size(res_MIQ.S, 1))};
        %res_GO = run_GO(fp, F0, V0, DEGREE);
        %res_TC_GO = TCODS(mesh, res_GO.S, F0, THETA0, DEGREE, VERBOSE);
        
        opt.cot = true;
        opt.graph = false;
        
        opt.gpu = true;
        opt.chol = false;
        opt.eigs = false;
        res_iter_gpu = run_lattice_iterative(mesh, ...
            F0, ...
            THETA0, ...
            DEGREE, ...
            N_ITER, ...
            opt);
        res_iter_gpu.title = {'Iterative (gpu)', ...
            sprintf('$E_{MIQ} = %g$', res_iter_gpu.E), ...
            sprintf('Elapsed time: %g', res_iter_gpu.elapsed_iter), ...
            sprintf('n sing: %d', size(res_iter_gpu.S, 1))};
        
%         opt.gpu = true;
%         opt.chol = false;
%         opt.eigs = true;
%         opt.n_eig = 100;
%         res_eigs = run_lattice_iterative(mesh, ...
%             F0, ...
%             THETA0, ...
%             DEGREE, ...
%             N_ITER, ...
%             opt);
%         res_eigs.title = {'Iterative (eigs)', ...
%             sprintf('$E_{MIQ} = %g$', res_eigs.E), ...
%             sprintf('Elapsed time: %g', res_eigs.elapsed_iter), ...
%             sprintf('n sing: %d', size(res_eigs.S, 1))};
        
        n_faces(k) = mesh.nF;
        
        T_MIQ(k) = res_MIQ.elapsed;
        T_iter_gpu(k) = res_iter_gpu.elapsed_iter;
        %T_eigs(k) = res_eigs.elapsed_iter;
        
        E_MIQ(k) = res_MIQ.E;
        E_iter_gpu(k) = res_iter_gpu.E;
        %E_eigs(k) = res_eigs.E;
        
        n_sing_MIQ(k) = size(res_MIQ.S, 1);
        n_sing_iter(k) = size(res_iter_gpu.S, 1);
        %n_sing_eigs(k) = size(res_eigs.S, 1);
        
        meshes{k} = mesh;
        results_MIQ{k} = res_MIQ;
        results_iter{k} = res_iter_gpu;
        %results_eigs{k} = res_eigs;
        
        k = k+1; 
    end
    
    figure(fig1);
    plot(n_faces, T_MIQ, '--xr', ...
        n_faces, T_iter_gpu, '--xb')
    legend('MIQ', ...
        'Iterative (gpu, nochol)')
    ylabel('Time (s)');
    xlabel('Number of faces');
    title('Time (s)')

    figure(fig2);
    plot(n_faces, E_MIQ, '--xr', ...
        n_faces, E_iter_gpu, '--xb')
    legend('MIQ', ...
        'Iterative (gpu, nochol)')
    ylabel('E');
    xlabel('Number of faces');
    title('Energy')

    figure(fig3);
    plot(n_faces, n_sing_MIQ, '--xr', ...
        n_faces, n_sing_iter, '--xb')
    legend('MIQ', ...
        'Iterative (gpu, nochol)')
    ylabel('Number of singularities');
    xlabel('Number of faces');
    title('Number of singularities')
    
    mysave(fig1, sprintf('%s_timing', out_filename));
    mysave(fig2, sprintf('%s_energy', out_filename));
    mysave(fig3, sprintf('%s_n_singularities', out_filename));
    
    MeshVis.wfigs([out_filename, '_cross_fields'], ...
                    {meshes{:}, meshes{:}}, ...
                    'OutFolder', OUT_FOLDER, ...
                    'Nrosy', {results_MIQ{:}, results_iter{:}}, ...
                    'ConstrainedFaces', F0, ...
                    'ConstraintVectors', V0, ...
                    'Save', SAVE, ...
                    'Rows', 2);
    
    
function timing_plot(meshes)
    global LOG;   
    global DEGREE;
    
    N_ITER = 1000;
    F0 = [1];
    THETA0 = [0];
    
    n_meshes = length(meshes);
    n_faces = zeros(n_meshes, 1);
    T_MIQ = zeros(n_meshes, 1);
    T_iter_gpu_chol = zeros(n_meshes, 1);
    T_iter_gpu_nochol = zeros(n_meshes, 1);
    T_iter_nogpu_chol = zeros(n_meshes, 1);
    T_iter_nogpu_nochol = zeros(n_meshes, 1);
    
    E_MIQ = zeros(n_meshes, 1);
    E_iter_gpu_chol = zeros(n_meshes, 1);
    E_iter_gpu_nochol = zeros(n_meshes, 1);
    E_iter_nogpu_chol = zeros(n_meshes, 1);
    E_iter_nogpu_nochol = zeros(n_meshes, 1);
    
    k = 1;
    for fname = meshes
        log_and_print(LOG, 'Loading %s...\r\n', fname{:});
        
        p = find_data_folder();
        fp = fullfile(p, fname{:});

        mesh = Mesh();
        mesh.loadTM(fp);   
        V0 = mesh.V(mesh.F(F0, 2), :) - mesh.V(mesh.F(F0, 1), :);   

        res_MIQ = run_MIQ(fp, F0, V0, DEGREE);
        %res_GO = run_GO(fp, F0, V0, DEGREE);
        %res_TC_GO = TCODS(mesh, res_GO.S, F0, THETA0, DEGREE, VERBOSE);
        
        opt.cot = true;
        opt.graph = false;
        
        opt.gpu = true;
        opt.chol = true;
        res_iter_gpu_chol = run_lattice_iterative(mesh, ...
           F0, ...
           THETA0, ...
           DEGREE, ...
           N_ITER, ...
           opt);
 
        opt.gpu = true;
        opt.chol = false;
        res_iter_gpu_nochol = run_lattice_iterative(mesh, ...
            F0, ...
            THETA0, ...
            DEGREE, ...
            N_ITER, ...
            opt);
        
        opt.gpu = false;
        opt.chol = true;
        res_iter_nogpu_chol = run_lattice_iterative(mesh, ...
           F0, ...
           THETA0, ...
           DEGREE, ...
           N_ITER, ...
           opt);
        
        opt.gpu = false;
        opt.chol = false;
        res_iter_nogpu_nochol = run_lattice_iterative(mesh, ...
            F0, ...
            THETA0, ...
            DEGREE, ...
            N_ITER, ...
            opt);
        
        n_faces(k) = mesh.nF;
        
        T_MIQ(k) = res_MIQ.elapsed;
        T_iter_gpu_chol(k) = res_iter_gpu_chol.elapsed;
        T_iter_gpu_nochol(k) = res_iter_gpu_nochol.elapsed;
        T_iter_nogpu_chol(k) = res_iter_nogpu_chol.elapsed;
        T_iter_nogpu_nochol(k) = res_iter_nogpu_nochol.elapsed;
        
        E_MIQ(k) = res_MIQ.E;
        E_iter_gpu_chol(k) = res_iter_gpu_chol.E;
        E_iter_gpu_nochol(k) = res_iter_gpu_nochol.E;
        E_iter_nogpu_chol(k) = res_iter_nogpu_chol.E;
        E_iter_nogpu_nochol(k) = res_iter_nogpu_nochol.E;
        
        k = k+1;

        %res_TC_GO.title = {'TC (w/ GO singularities)', ...
        %    sprintf('$E_{MIQ} = %g$', res_TC_GO.E)};

        %results = {res_MIQ, res_iterative, res_GO, res_TC_GO};
        %results = {res_iterative};
%         results = {res_MIQ, res_iterative};
% 
%         MeshVis.wfigs(fname{:}, ...
%                     mesh, ...
%                     'OutFolder', OUT_FOLDER, ...
%                     'Nrosy', results, ...
%                     'ConstrainedFaces', F0, ...
%                     'ConstraintVectors', V0, ...
%                     'Save', SAVE);
    end
    
    figure;
    plot(n_faces, T_MIQ, '--xr', ...
        n_faces, T_iter_gpu_chol, '--xb', ...
        n_faces, T_iter_gpu_nochol, '--xg', ...
        n_faces, T_iter_nogpu_chol, '--xm', ...
        n_faces, T_iter_nogpu_nochol, '--xk')
    legend('MIQ', ...
        'Iterative (gpu, chol)', ...
        'Iterative (gpu, nochol)', ...
        'Iterative (nogpu, chol)', ...
        'Iterative (nogpu, nochol)');
    ylabel('Time (s)');
    xlabel('Number of faces');
    
    figure;
    plot(n_faces, E_MIQ, '--xr', ...
        n_faces, E_iter_gpu_chol, '--xb', ...
        n_faces, E_iter_gpu_nochol, '--xg', ...
        n_faces, E_iter_nogpu_chol, '--xm', ...
        n_faces, E_iter_nogpu_nochol, '--xk')
    legend('MIQ', ...
        'Iterative (gpu, chol)', ...
        'Iterative (gpu, nochol)', ...
        'Iterative (nogpu, chol)', ...
        'Iterative (nogpu, nochol)');
    ylabel('E');
    xlabel('Number of faces');

function lap_pseudo_inverse_timings(mesh_filenames)
    global LOG;
    global EPS;
    global DEGREE;
    
    N_ITER = 1000;
    
    data_folder = find_data_folder();
    
    F0 = [1];
    THETA0 = [0];

    r = 1;
    for mname = mesh_filenames
        gd = gpuDevice(1);
        fp = fullfile(data_folder, mname{:});
        disp(fp)
    
        mesh = Mesh();
        mesh.loadTM(fp);  
        nV = mesh.nV;

        L = lap_cot(mesh);

        %tic; Lp_vanilla    = pinv(full(L))                           ; elapsed_vanilla(r) = toc;
        %tic; Lp_spqr_pinv  = spqr_pinv(L, speye(size(L)))            ; elapsed_spqr_pinv(r) = toc;
        %tic; Lp_qrginv     = qrginv(L)                               ; elapsed_qrginv(r) = toc;
        tic; Lp_chol        = invChol_mex(L + 1/nV) - 1/nV            ; elapsed_chol(r) = toc;
        tic; Lp_gpu_sing   = inv(gpuArray(single(L + 1/nV))) - 1/nV  ; wait(gd); elapsed_gpu_sing(r) = toc;
        %tic; Lp_gpu_dbl    = inv(gpuArray((L + 1/nV))) - 1/nV        ; wait(gd); elapsed_gpu_dbl(r) = toc;
        tic; Lp_d_and_c    = pinv_lap_divide_and_conq(L)             ; wait(gd); elapsed_d_and_c(r) = toc;
        %elapsed_gpu_sing(r) = gputimeit(@() (inv(gpuArray(single(L+1/nV)))-1/nV));
        %elapsed_gpu_dbl(r)  = gputimeit(@() (inv(gpuArray(L+1/nV)) - 1/nV));
        %elapsed_d_and_c(r)  = gputimeit(@() pinv_lap_divide_and_conq(L));
        %tic; Lp_spqr_pinv_gpu = spqr_pinv(gpuArray(full(L)), eye(size(L))); elapsed_spqr_pinv_gpu = toc;

        % Iterative improvement loop
        gd = gpuDevice(1);
        Lp = (inv(gpuArray(single(L+1/nV))) - 1/nV);
        tic
        Ad = get_gaussian_curvature(mesh);
        k0 = (2/pi)*Ad;
        c = round(sum(k0));

        k = zeros(nV, 1);
        inds = randperm(nV, c);
        k(inds) = 1;

        % Resistance distance matrix
        dlp = diag(Lp); 
        R = bsxfun(@plus,dlp,dlp') - 2*Lp; 
        b = gpuArray(2*Lp*(k-k0));

        for iter = 1:N_ITER
            D = bsxfun(@minus, b, b');

            Z = D + R;
            [m, i] = min(Z(:));

            if abs(m) > EPS
                break;
            end

            [i, j] = ind2sub(size(Z), i);
            k(i) = k(i) + 1;
            k(j) = k(j) - 1;

            b = b + 2*Lp(:,i) - 2*Lp(:,j);

            clear Z;
        end
        wait(gd); elapsed_iter(r) = toc;
        
        opt.cot = true;
        opt.graph = false;
        
        opt.gpu = true;
        opt.chol = false;
        opt.eigs = false;
        res_iter_gpu = run_lattice_iterative(mesh, ...
            F0, ...
            THETA0, ...
            DEGREE, ...
            N_ITER, ...
            opt);
        elapsed_lat_iter(r) = res_iter_gpu.elapsed_iter;
        
        V0 = mesh.V(mesh.F(F0, 2), :) - mesh.V(mesh.F(F0, 1), :);   
        res_MIQ = run_MIQ(fp, F0, V0, DEGREE);
        elapsed_MIQ(r) = res_MIQ.elapsed;
        
        % Approximate inverse using eigendecomposition
%         tic
%         n_eig = 100;
%         [V, D] = eigs(L, n_eig, 'sm');
%         Dp = diag(D);
%         Dp(Dp < EPS) = inf;
%         Dp = 1 ./ Dp;
%         Dp = diag(Dp);
%         Lp_eigs = V*Dp*V';
%         elapsed_eigs(k) = toc;
        
        n_faces(r) = mesh.nF;

        %log_and_print(LOG, 'pinv(full(L))  : %g\r\n', elapsed_vanilla(k));
        %log_and_print(LOG, 'spqr_pinv      : %g\r\n', elapsed_spqr_pinv(k));
        %log_and_print(LOG, 'qrginv(L)      : %g\r\n', elapsed_qrginv(k));
        log_and_print(LOG, 'chol mex       : %g\r\n', elapsed_chol(r));
        log_and_print(LOG, 'gpu, single    : %g\r\n', elapsed_gpu_sing(r));
        %log_and_print(LOG, 'gpu, double    : %g\r\n', elapsed_gpu_dbl(r));
        log_and_print(LOG, 'Div and conq   : %g\r\n', elapsed_d_and_c(r));
        log_and_print(LOG, 'elapsed iter   : %g\r\n', elapsed_iter(r));
        log_and_print(LOG, 'lattice iter   : %g\r\n', elapsed_lat_iter(r));
        log_and_print(LOG, 'MIQ            : %g\r\n', elapsed_MIQ(r));
        %log_and_print(LOG, 'eigs           : %g\r\n', elapsed_eigs(k));
        
        r = r + 1;
    end
    
    fig = figure();    
    %plot(...
        %n_faces, elapsed_vanilla, '--xr', ...
        %n_faces, elapsed_spqr_pinv, '--xb', ...
        %n_faces, elapsed_qrginv, '--xm', ...
     plot(   n_faces, elapsed_chol, '--xg', ...
        n_faces, elapsed_gpu_sing, '--xy', ...
        n_faces, elapsed_d_and_c, '--xk', ...
        n_faces, elapsed_iter, '--r', ...
        n_faces, elapsed_lat_iter, '--m', ...
        n_faces, elapsed_MIQ', '--b')
        %n_faces, elapsed_gpu_dbl, '--xc', ...
    %legend(...
        %'pinv', ...
        %'spqr pinv', ...
        %'qrginv', ...
    legend(    'inv chol mex', ...
        'inv gpu (single)', ...
        'inv d&c partition step', ...
        'lattice (loop only)', ...
        'lattice (total)', ...
        'MIQ')
        %'inv gpu (double)', ...
    ylabel('Time (s)');
    xlabel('Number of faces');
    title('Laplacian inverse timing')
    
    mysave(fig, 'lap_inv_timing');
    
function energy_and_time_plots(meshes)
    global LOG;   
    global DEGREE;
    
    N_ITER = 1000;
    N_REPEAT = 10;
    F0 = [1];
    THETA0 = [0];
    
    n_meshes = length(meshes);
    n_faces = zeros(n_meshes, 1);
    T_MIQ = zeros(n_meshes, N_REPEAT);
    T_iter_gpu_chol = zeros(n_meshes, N_REPEAT);
    T_iter_gpu_nochol = zeros(n_meshes, N_REPEAT);
    T_iter_nogpu_chol = zeros(n_meshes, N_REPEAT);
    T_iter_nogpu_nochol = zeros(n_meshes, N_REPEAT);
    
    E_MIQ = zeros(n_meshes, N_REPEAT);
    E_iter_gpu_chol = zeros(n_meshes, N_REPEAT);
    E_iter_gpu_nochol = zeros(n_meshes, N_REPEAT);
    E_iter_nogpu_chol = zeros(n_meshes, N_REPEAT);
    E_iter_nogpu_nochol = zeros(n_meshes, N_REPEAT);
    
    k = 1;
    for fname = meshes
        for r = 1:N_REPEAT
            log_and_print(LOG, 'Loading %s...\r\n', fname{:});

            p = find_data_folder();
            fp = fullfile(p, fname{:});

            mesh = Mesh();
            mesh.loadTM(fp);   
            V0 = mesh.V(mesh.F(F0, 2), :) - mesh.V(mesh.F(F0, 1), :);   

            res_MIQ = run_MIQ(fp, F0, V0, DEGREE);
            %res_GO = run_GO(fp, F0, V0, DEGREE);
            %res_TC_GO = TCODS(mesh, res_GO.S, F0, THETA0, DEGREE, VERBOSE);

            opt.cot = true;
            opt.graph = false;
            
            opt.gpu = true;
            opt.chol = false;
            res_iter_gpu_nochol = run_lattice_iterative(mesh, ...
                F0, ...
                THETA0, ...
                DEGREE, ...
                N_ITER, ...
                opt);

            opt.gpu = true;
            opt.chol = true;
            res_iter_gpu_chol = run_lattice_iterative(mesh, ...
               F0, ...
               THETA0, ...
               DEGREE, ...
               N_ITER, ...
               opt);

            opt.gpu = false;
            opt.chol = true;
            res_iter_nogpu_chol = run_lattice_iterative(mesh, ...
               F0, ...
               THETA0, ...
               DEGREE, ...
               N_ITER, ...
               opt);

            opt.gpu = false;
            opt.chol = false;
            res_iter_nogpu_nochol = run_lattice_iterative(mesh, ...
                F0, ...
                THETA0, ...
                DEGREE, ...
                N_ITER, ...
                opt);

            n_faces(k) = mesh.nF;

            T_MIQ(k, r) = res_MIQ.elapsed;
            T_iter_gpu_chol(k, r) = res_iter_gpu_chol.elapsed_iter;
            T_iter_gpu_nochol(k, r) = res_iter_gpu_nochol.elapsed_iter;
            T_iter_nogpu_chol(k, r) = res_iter_nogpu_chol.elapsed_iter;
            T_iter_nogpu_nochol(k, r) = res_iter_nogpu_nochol.elapsed_iter;

            E_MIQ(k, r) = res_MIQ.E;
            E_iter_gpu_chol(k, r) = res_iter_gpu_chol.E;
            E_iter_gpu_nochol(k, r) = res_iter_gpu_nochol.E;
            E_iter_nogpu_chol(k, r) = res_iter_nogpu_chol.E;
            E_iter_nogpu_nochol(k, r) = res_iter_nogpu_nochol.E;

            %res_TC_GO.title = {'TC (w/ GO singularities)', ...
            %    sprintf('$E_{MIQ} = %g$', res_TC_GO.E)};

            %results = {res_MIQ, res_iterative, res_GO, res_TC_GO};
            %results = {res_iterative};
    %         results = {res_MIQ, res_iterative};
    % 
    %         MeshVis.wfigs(fname{:}, ...
    %                     mesh, ...
    %                     'OutFolder', OUT_FOLDER, ...
    %                     'Nrosy', results, ...
    %                     'ConstrainedFaces', F0, ...
    %                     'ConstraintVectors', V0, ...
    %                     'Save', SAVE);
        end
        
        k = k+1;
    end
    
    figure;
    hold on
    plot(n_faces, mean(T_MIQ, 2), '--xr', ...
        n_faces, mean(T_iter_gpu_chol, 2), '--xb', ...
        n_faces, mean(T_iter_gpu_nochol, 2), '--xg', ...
        n_faces, mean(T_iter_nogpu_chol, 2), '--xm', ...
        n_faces, mean(T_iter_nogpu_nochol, 2), '--xk')
    errorbar(n_faces, mean(T_MIQ, 2), std(T_MIQ, 0, 2), '--xr')
    errorbar(n_faces, mean(T_iter_gpu_chol, 2), std(T_iter_gpu_chol, 0, 2), '--xb')
    errorbar(n_faces, mean(T_iter_gpu_nochol, 2), std(T_iter_gpu_nochol, 0, 2), '--xg')
    errorbar(n_faces, mean(T_iter_nogpu_chol, 2), std(T_iter_nogpu_chol, 0, 2), '--xm')
    errorbar(n_faces, mean(T_iter_nogpu_nochol, 2), std(T_iter_nogpu_nochol, 0, 2), '--xk')
    hold off
    legend('MIQ', ...
        'Iterative (gpu, chol)', ...
        'Iterative (gpu, nochol)', ...
        'Iterative (nogpu, chol)', ...
        'Iterative (nogpu, nochol)');
    ylabel('Time (s)');
    xlabel('Number of faces');
    
    figure;
    hold on
    plot(n_faces, mean(E_MIQ, 2), '--xr', ...
        n_faces, mean(E_iter_gpu_chol, 2), '--xb', ...
        n_faces, mean(E_iter_gpu_nochol, 2), '--xg', ...
        n_faces, mean(E_iter_nogpu_chol, 2), '--xm', ...
        n_faces, mean(E_iter_nogpu_nochol, 2), '--xk')
    errorbar(n_faces, mean(E_MIQ, 2), std(E_MIQ, 0, 2), '--xr')
    errorbar(n_faces, mean(E_iter_gpu_chol, 2), std(E_iter_gpu_chol, 0, 2), '--xb')
    errorbar(n_faces, mean(E_iter_gpu_nochol, 2), std(E_iter_gpu_nochol, 0, 2), '--xg')
    errorbar(n_faces, mean(E_iter_nogpu_chol, 2), std(E_iter_nogpu_chol, 0, 2), '--xm')
    errorbar(n_faces, mean(E_iter_nogpu_nochol, 2), std(E_iter_nogpu_nochol, 0, 2), '--xk')
    hold off
    legend('MIQ', ...
        'Iterative (gpu, chol)', ...
        'Iterative (gpu, nochol)', ...
        'Iterative (nogpu, chol)', ...
        'Iterative (nogpu, nochol)');
    ylabel('Energy');
    xlabel('Number of faces');

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
    
    res.elapsed = str2double(elapsed{:}{:});
    
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
  
function gpu_memory_tests(fp)
    global LOG;
    global DEGREE;
    gd = gpuDevice(1);
    
    F0 = [1];
    THETA0 = [0];
    N_ITER = 1000;
    
    mesh = Mesh();
    mesh.loadTM('../data/bunny_57k_faces.off');
    
    %V0 = mesh.V(mesh.F(F0, 2), :) - mesh.V(mesh.F(F0, 1), :);   
    %res_MIQ = run_MIQ('../data/bunny_57k_faces.off', F0, V0, DEGREE);
    %elapsed_MIQ = res_MIQ.elapsed;
    
    opt.cot = true;
    opt.graph = false;
    opt.gpu = false;
    opt.only_Lp_on_gpu = true;
    res_iter_gpu = run_lattice_iterative(mesh, ...
        F0, ...
        THETA0, ...
        DEGREE, ...
        N_ITER, ...
        opt);
    elapsed_lat_iter(r) = res_iter_gpu.elapsed_iter;

    mesh = Mesh();
    mesh.loadTM(fp);
    
    nF = mesh.nF;
    nV = mesh.nV;
    Ad = get_gaussian_curvature(mesh);
    k0 = (2/pi)*Ad;
    c = round(sum(k0));
    k = zeros(nV, 1);
    inds = randperm(nV, c);
    k(inds) = 1;
    
    log_and_print(LOG, 'gpu_memory_tests...\n');
    log_and_print(LOG, 'nF : %d\n', nF);
    log_and_print(LOG, 'nV : %d\n', nV);
    wait(gd); avail_mem = gd.AvailableMemory / 1024^2;
    log_and_print(LOG, 'Available memory: %g\n\n', avail_mem);
    
    log_and_print(LOG, 'Creating L+1/nV...\n');
    L = gpuArray(single(lap_cot(mesh)-1/nV));
    wait(gd); avail_mem = gd.AvailableMemory / 1024^2;
    log_and_print(LOG, 'Available memory: %g\n\n', avail_mem);
    
    log_and_print(LOG, 'Inverting L...\n');
    Lp = inv(L) - 1/nV;
    wait(gd); avail_mem = gd.AvailableMemory / 1024^2;
    log_and_print(LOG, 'Available memory: %g\n\n', avail_mem);
    
    log_and_print(LOG, 'Clearing L...\n');
    clear L;
    wait(gd); avail_mem = gd.AvailableMemory / 1024^2;
    log_and_print(LOG, 'Available memory: %g\n\n', avail_mem);
    
    log_and_print(LOG, 'Creating R...\n');
    dlp = diag(Lp); 
    R = bsxfun(@plus,dlp,dlp') - 2*Lp;
    wait(gd); avail_mem = gd.AvailableMemory / 1024^2;
    log_and_print(LOG, 'Available memory: %g\n\n', avail_mem);
    
    log_and_print(LOG, 'Clearing dlp...\n');
    clear dlp;
    wait(gd); avail_mem = gd.AvailableMemory / 1024^2;
    log_and_print(LOG, 'Available memory: %g\n\n', avail_mem);
    
    log_and_print(LOG, 'Creating b...\n');
    b = gpuArray(2*Lp*(k-k0));
    wait(gd); avail_mem = gd.AvailableMemory / 1024^2;
    log_and_print(LOG, 'Available memory: %g\n\n', avail_mem);
    
    log_and_print(LOG, 'Creating D...\n');
    D = bsxfun(@minus, b, b');
    wait(gd); avail_mem = gd.AvailableMemory / 1024^2;
    log_and_print(LOG, 'Available memory: %g\n\n', avail_mem);
    
    log_and_print(LOG, 'Creating Z...\n');
    Z = D + R;
    [m, i] = min(Z(:));
    wait(gd); avail_mem = gd.AvailableMemory / 1024^2;
    log_and_print(LOG, 'Available memory: %g\n\n', avail_mem);
    
function test_CUDA(fp)
    global LOG;
    global DEGREE;
    gd = gpuDevice(1);
    
    F0 = [1];
    THETA0 = [0];
    N_ITER = 1000;
    
    mesh = Mesh();
    mesh.loadTM(fp);
    
    opt.cot = true;
    opt.graph = false;
    opt.gpu = false;
    opt.CUDAKernel = true;
    res_CUDA = run_lattice_iterative(mesh, ...
        F0, ...
        THETA0, ...
        DEGREE, ...
        N_ITER, ...
        opt);
    
    
function res = run_lattice_iterative(mesh, f0, theta0, degree, n_iter, opt)
    % Iterative {+1,-1} improvement heuristic.
    % See notes/2017-08-31 Iterative impr[3852].pdf
    %
    % Input:
    %   mesh - the mesh on which to calculate the direction field
    %   f0 - constrained faces
    %   theta0 - constraint angles
    %   degree - field degree
    %   n_iter - number of iterations
    %   opt - struct with the following boolean fields:
    %           graph - use graph Laplacian
    %           cot - use cot Laplacian
    %           gpu - use GPU
    %           chol - use cholesky decomposition
    %           eigs - use eigen decomposition
    %
    % Output:
    %   res - struct with the following fields:
    %           elapsed_iter - our heuristic elapsed time (in seconds)
    %           elapsed_TC - TC elapsed time
    %           elapsed - total elapsed time
    %           see TCODS for the rest.
    
    if ~isfield(opt, 'plot') , opt.plot = false  ; end
    if ~isfield(opt, 'gpu')  , opt.gpu = false; opt.only_Lp_on_gpu=false ; end
    if ~isfield(opt, 'graph'), opt.graph = false ; opt.cot = true; end
    if ~isfield(opt, 'cot')  , opt.graph = true  ; opt.cot = false; end
    if ~isfield(opt, 'chol') , opt.chol = false  ; end
    if ~isfield(opt, 'eigs') , opt.eigs = false  ; end
    if ~isfield(opt, 'only_Lp_on_gpu') , opt.only_Lp_on_gpu = false ; end
    if ~isfield(opt, 'CUDAKernel') , opt.CUDAKernal = false ; end
    if ~isfield(opt, 'CUDAKernel_reduce_cols'), opt.CUDAKernel_reduce_cols = false ; end
    if ~isfield(opt, 'TCODS'), opt.TCODS = false; end
       
    if opt.gpu || opt.only_Lp_on_gpu || opt.CUDAKernel
        gpuDevice(1);
    end
    
    tic
    global EPS;
    global LOG;
    
    nV = mesh.nV;
    
    log_and_print(LOG, 'Creating L...\r\n');
    if opt.graph
        [d0, ~] = get_exterior_derivatives(mesh);
        L = d0'*d0;
    elseif opt.cot
        %L = lap(mesh);
        L = lap_cot(mesh);
    end
    
    elapsed_setup = toc();
    
    %xi = 2 - 2*mesh.genus;
     
    % faster way to get the pseudo inverse   
    log_and_print(LOG, 'Inverting L...\r\n');
    tic
    if opt.chol
        Lp = invChol_mex(L + 1/nV) - 1/nV;
        if opt.gpu
            Lp = gpuArray(single(Lp));
        end
    elseif opt.eigs
        [V, Dp] = eigs(L, opt.n_eig, 'sm');
        Dp = diag(Dp);
        Dp(Dp < EPS) = inf;
        Dp = diag(1 ./ Dp);
        Lp = V*Dp*V';
        Lp = gpuArray(full(single(Lp)));
    elseif opt.gpu
        Lp = inv(gpuArray(single(L + 1/nV))) - 1/nV;
    elseif opt.only_Lp_on_gpu
        Lp = gpuArray(invChol_mex(single(L+1/nV)) - 1/nV);
    else
        Lp = inv(L + 1/nV) - 1/nV;
    end
    elapsed_inv = toc();
    
    % Rather disappointing
    %Lp = qrginv(L);

    tic
    Ad = get_gaussian_curvature(mesh);
    k0 = (2/pi)*Ad;
    c = round(sum(k0)); % = 4*xi
    
    k = zeros(nV, 1);
    inds = randperm(nV, c);
    k(inds) = 1;
    
    % Resistance distance matrix
    log_and_print(LOG, 'Calculating resistance matrix...\r\n');
    
    if opt.gpu
        % Not clear which one is better
        %Ldii = repmat(diag(Lp), 1, nV);
        %R = gpuArray(Ldii + Ldii' - 2*Lp);
        %tic
        dlp = diag(Lp); 
        R = bsxfun(@plus,dlp,dlp') - 2*Lp; 
        %toc
        b = gpuArray(2*Lp*(k-k0)); % runs out of memory with backslash (on large meshes)
    elseif opt.only_Lp_on_gpu
        error('Doesn''t work')
        dlp = diag(Lp); 
        %R = gather(bsxfun(@plus,dlp,dlp') - 2*Lp); 
        b = 2*(Lp*(k-k0)); % runs out of memory with backslash (on large meshes)
    elseif opt.CUDAKernel
        % Setup the kernel
        kernel = parallel.gpu.CUDAKernel('lattice_inner_loop.ptx', 'lattice_inner_loop.cu', 'lattice_inner_loop');
        num_elem = nV^2;
        kernel.ThreadBlockSize = [kernel.MaxThreadsPerBlock, 1, 1];
        kernel.GridSize = [ceil(num_elem / kernel.MaxThreadsPerBlock), 1];
        
        dlp = diag(Lp);
        b = gpuArray(2*Lp*(k-k0));
        out = zeros(size(Lp), 'like', Lp);
    elseif opt.CUDAKernel_reduce_cols
        % Setup the kernel
        kernel = parallel.gpu.CUDAKernel('lattice_inner_loop.ptx', 'lattice_inner_loop.cu', 'reduce_cols');
        num_elem = nV^2;
        kernel.ThreadBlockSize = [1, 1, 1];
        kernel.GridSize = [1, nV];
        
        dlp = diag(Lp);
        b = gpuArray(2*Lp*(k-k0));
        out_min = zeros(1, nV, 'like', Lp);
        out_r = zeros(1, nV, 'like', Lp);
    else
        %Ldii = repmat(diag(Lp), 1, nV);
        %R = gather(Ldii + Ldii' - 2*Lp);
        %b = 2*Lp*(k-k0);
        dlp = diag(Lp); 
        R = bsxfun(@plus,dlp,dlp') - 2*Lp; 
        b = 2*(Lp \ (k-k0));
    end
    elapsed_setup = elapsed_setup + toc();
    
    log_and_print(LOG, 'Iterative improvement...\r\n');
    tic
    for iter = 1:n_iter
        if opt.only_Lp_on_gpu
            error('Doesn''t work')
            [m, i] = min(bsxfun(@plus, dlp, dlp') - 2*Lp + bsxfun(@minus, b, b'));
            [m, j] = min(m);
            i = i(j);
        elseif opt.CUDAKernel
            out = feval(kernel, out, Lp, dlp, b, nV, nV, num_elem);
            [m, i] = min(Z(:));
            [i, j] = ind2sub(size(Z), i);
        elseif opt.CUDAKernel_reduce_cols
            [out_min, out_r] = feval(kernel, out_min, out_r, Lp, dlp, b, nV, nV, num_elem);
            [m, j] = min(out_min);
            i = out_r(j);
        else
            D = bsxfun(@minus, b, b');
            Z = D + R;
            [m, i] = min(Z(:));
            [i, j] = ind2sub(size(Z), i);
        end
        
        if abs(m) < EPS
            break;
        end

        k(i) = k(i) + 1;
        k(j) = k(j) - 1;
        
        b = b + 2*Lp(:,i) - 2*Lp(:,j);
        
        if opt.plot
            E_lat_hist(iter) = gather((k-k0)'*Lp*(k-k0));
            m_hist(iter) = gather(m);
            min_D_hist(iter) = gather(min(D(:)));
            min_R_hist(iter) = gather(min(R(:)));
        end
        
        if ~opt.only_Lp_on_gpu
            clear Z;
        end
    end
    wait(gpuDevice); elapsed_iter = toc;
    
    if opt.plot
        figure
        xx = 1:length(E_lat_hist);
        plot(xx, E_lat_hist, '-r', ...
            xx, m_hist, '-b', ...
            xx, min_D_hist, '-m', ...
            xx, min_R_hist, '-g')
        legend('E lat', ...
            'm', ...
            'min(D)', ...
            'min(R)')
    end
    
%     tic
    if opt.TCODS
        inds = find(abs(k) > EPS);
        S = [inds, k(inds) / degree];
        res = TCODS(mesh, S, f0, theta0, degree, false);
    end
     
     res.S = S;
     res.k = k;
     res.elapsed_setup = elapsed_setup;
     res.elapsed_inv = elapsed_inv;
     res.elapsed_iter = elapsed_iter;
     res.elapsed_total = elapsed_setup + elapsed_inv + elapsed_iter;
     
%     res.elapsed_TC = toc;
%     res.elapsed_iter = elapsed_iter;
%     
%     log_and_print(LOG, 'E = %g\r\n', res.E);
%     log_and_print(LOG, 'Done.\r\n');
    
%     res.elapsed = res.elapsed_iter + res.elapsed_TC;
    res.title = {'Iterative', ...
        sprintf('Setup time: %g', res.elapsed_setup), ...
        sprintf('Inv time: %g', res.elapsed_inv), ...
        sprintf('Iter time: %g', res.elapsed_iter), ...
        sprintf('Total time: %g', res.elapsed_total)};
    if opt.TCODS
        res.title = {res.title{:}, sprintf('$E_{MIQ} = %g$', res.E)};
    end

    
function ginv = ginv(X)
    %Returns the Moore-Penrose inverse of the argument
    if isempty(X)
        quick return
        ginv = zeros(size(X'), class(X));
        return
    end
    
    [n,m] = size(X);
    if sprank(X) < min(n,m);
        error('matrix must be of full rank');
    elseif n > m,
        C = X' * X ;
        ginv = C \ X';
    else
        C = X * X';
        G = C \ X;
        ginv = G';
    end
    
function qrginv = qrginv(B) 
    [N,M] = size(B); 
    [Q,R,P] = spqr(B);
    r=sum(any(abs(R)>1e-5,2)); 
    R1 = R(1:r,:); 
    R2 = ginv(R1); 
    R3 = [R2 zeros(M,N-r)]; 
    A = P*R3*Q'; 
    qrginv = A;

function Lp = pinv_lap_divide_and_conq(L)
    assert(size(L, 1) == size(L, 2))
    n = size(L, 1);

    % Partition the graph
    [V, ~] = eigs(L, 2, 'SA');
    %inds = V(:, 2) > 0;
    %L1 = gpuArray(full(L(inds, inds)));
    %L2 = gpuArray(full(L(~inds, ~inds)));
    
    % Invert the partitions
    %Lp1 = inv(L1 + 1/n) - 1/n;
    %Lp2 = inv(L2 + 1/n) - 1/n;
    
    Lp = [];
    
    
function subdivide_and_save(fp, n_sbd)
    mesh = Mesh();
    mesh.loadTM(fp);
    [pathstr, name, ext] = fileparts(fp);
    mesh.saveTM(fullfile(pathstr, sprintf('%s_sd0%s', name, ext)));
    for i = 1:n_sbd
        mesh = subdivide_mesh(mesh);
        new_fp = fullfile(pathstr, sprintf('%s_sd%d%s', name, i, ext));
        mesh.saveTM(new_fp);
    end
    
function new_mesh = subdivide_mesh(mesh)
    EV = mesh.EVAdj;
    nE = mesh.nE;
    nV = mesh.nV;
    nF = mesh.nF;
    V = mesh.V;
    F = mesh.F;
    V_new = zeros(nE, 3);
    
    VV_to_eids = sparse([], [], [], mesh.nV, mesh.nV, mesh.nE);
    for eid = 1:mesh.nE
        v1 = EV(eid, 1);
        v2 = EV(eid, 2);
        assert(v1 < v2);
        VV_to_eids(v1, v2) = eid;
        VV_to_eids(v2, v1) = eid;
    end
    
    for eid = 1:nE
        v1 = V(EV(eid, 1), :);
        v2 = V(EV(eid, 2), :);
        V_new(eid, :) = (v1 + v2) / 2;
    end
    
    F_new = zeros(4*nF, 3);
    for fid = 1:nF
        v1 = F(fid, 1);
        v2 = F(fid, 2);
        v3 = F(fid, 3);
        
        e12 = VV_to_eids(v1, v2);
        e23 = VV_to_eids(v2, v3);
        e31 = VV_to_eids(v3, v1);
        
        v12 = nV + e12;
        v23 = nV + e23;
        v31 = nV + e31;
        
        idx = 4*(fid-1);
        F_new(idx+1, :) = [v1 , v12, v31];
        F_new(idx+2, :) = [v12, v2 , v23];
        F_new(idx+3, :) = [v23, v3 , v31];
        F_new(idx+4, :) = [v31, v12, v23];
    end
    
    V_new = [V; V_new];
    new_mesh = Mesh(V_new, F_new);
    
function mysave(fig, filename)
    global RESOLUTION;
    global ASPECT_RATIO;
    global OUT_FOLDER;
    
    %pngname = sprintf('%s_%04d.png', filename, i);
    filename = fullfile(OUT_FOLDER, filename);
    pngname = sprintf('%s.png', filename);

    dpi = get(0, 'ScreenPixelsPerInch');
    in = RESOLUTION/dpi;

    %fig = gcf;
    figure(fig);
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 ASPECT_RATIO*in in];
    fig.PaperPositionMode = 'manual';
    print(pngname,'-dpng','-r0')
C:\miri\Dropbox (Technion CS)\miri\technion\trips\obergurgl17C:\miri\Dropbox (Technion CS)\miri\technion\trips\obergurgl17