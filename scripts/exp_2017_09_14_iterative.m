function [] = exp_2017_08_31_1145_iterative()

    global LOG;
    global EPS;
    global VERBOSE;
    global PLOT;
    global SAVE;
    global DEGREE;
    global ASPECT_RATIO;
    global RESOLUTION;
    VERBOSE = true;
    EPS = 1e-9;
    PLOT = true;
    SAVE = true;
    DEGREE = 4;
    ASPECT_RATIO = 1;
    RESOLUTION = 1024;    
    LOG = -1;
    
    F0 = [1];
    THETA0 = [0];
    V0 = [1, 0, 0];
    N_ITER = 1000;
    
    NUMBER_OF_SINGULARITIES.run = false;
    
    MIQ_VS_OTC_BUNNIES_LARGE.name = ...
        'MIQ_vs_OTC_inputmodels_9gb_gpu_matlab2017b';
    MIQ_VS_OTC_BUNNIES_LARGE.description = ...
        'MIQ vs OTC (cpu inv, gpu loop) vs OTC (cpu inv, cpu loop) vs OTC (block inv, gpu loop) on bunnies.';
    MIQ_VS_OTC_BUNNIES_LARGE.run = true;
    MIQ_VS_OTC_BUNNIES_LARGE.merge = true;
    MIQ_VS_OTC_BUNNIES_LARGE.plot = true;
    data_folder = '../data/inputmodels';
    ext = '.obj';
    
    opt_vanilla = struct('inv_gpu', false, ...
        'iter_gpu', false, ...
        'TCODS', true);
    opt_inv_cpu = struct('inv_gpu', false, ...
        'CUDAKernel_reduce_cols', true, ...
        'TCODS', true);
    opt_block_inv_gpu = struct('block_inv_gpu', true, ...
        'CUDAKernel_reduce_cols', true, ...
        'TCODS', true, ...
        'lb', 20000);
    opt_inv_gpu = struct('inv_gpu', true, ...
        'CUDAKernel_reduce_cols', true, ...
        'TCODS', true);
    experiments = {
        struct(...
            'func', @run_lattice_iter, ...
            'args', {{F0, THETA0, DEGREE, N_ITER, opt_vanilla}}, ...
            'description', 'OTC (cpu inv, cpu loop)', ...
            'rerun', false), ...
        struct(...
            'func', @run_lattice_iter, ...
            'args', {{F0, THETA0, DEGREE, N_ITER, opt_inv_cpu}}, ...
            'description', 'OTC (cpu inv, gpu loop)', ...
            'rerun', false), ...
        struct(...
            'func', @run_lattice_iter, ...
            'args', {{F0, THETA0, DEGREE, N_ITER, opt_block_inv_gpu}}, ...
            'description', 'OTC (block inv, gpu loop)', ...
            'rerun', true), ...
        struct(...
            'func', @run_lattice_iter, ...
            'args', {{F0, THETA0, DEGREE, N_ITER, opt_inv_gpu}}, ...
            'description', 'OTC (gpu inv, gpu loop)', ...
            'rerun', false), ...
        struct(...
            'func', @run_MIQ, ...
            'args', {{F0, V0, DEGREE}}, ...
            'description', 'MIQ', ...
            'rerun', false), ...
        };
    
    plots = {...
        struct(...
            'title', 'Energy', ...
            'name', 'energy', ...
            'get_x_field', @(x) x.nF, ...
            'get_y_field', @(x) x.E, ...
            'style', {{'--.'}}), ...
        struct(...
            'title', 'Elapsed time (total)', ...
            'name', 'elapsed_total', ...
            'get_x_field', @(x) x.nF, ...
            'get_y_field', @get_elapsed_total, ...
            'style', {{'--.'}}), ...
        };
           
    if MIQ_VS_OTC_BUNNIES_LARGE.run || MIQ_VS_OTC_BUNNIES_LARGE.plot
        experiment_name = MIQ_VS_OTC_BUNNIES_LARGE.name;
        experiment_description = MIQ_VS_OTC_BUNNIES_LARGE.description;
        merge = MIQ_VS_OTC_BUNNIES_LARGE.merge;
        
        if MIQ_VS_OTC_BUNNIES_LARGE.run
            filepaths = get_filepaths(data_folder, ext);
            filepaths = filter_meshes(filepaths, @(m) abs(m.genus)<1e-10);
            %filepaths = {'../data/bunnies/bunny_57k_faces.off'};
            %filepaths = {'../data/bunnies/bunny_res1.off'};
            %filepaths = get_filepaths('../data/bunnies_small', ext);

            [results, out_folder] = run_experiments(...
                filepaths, ...
                experiments, ...
                experiment_name, ...
                experiment_description);
            
            %sendmail('nahumf@cs.technion.ac.il', sprintf('Done: %s', experiment_name));
            
            if merge
                out_folder = fullfile('..', 'results', 'experiments', experiment_name);
                if ~exist(out_folder, 'dir')
                    mkdir(out_folder);
                end
                old_results_fp = fullfile(out_folder, 'results.mat');
                if exist(old_results_fp, 'file')
                    old_results = load(old_results_fp);
                    results = merge_results(old_results.results, results, experiments, out_folder);
                else
                    save(old_results_fp, 'results');
                end
            end
        else
            out_folder = fullfile('..', 'results', 'experiments', experiment_name);
            load(fullfile(out_folder, 'results.mat'))
        end
        
        if MIQ_VS_OTC_BUNNIES_LARGE.plot
            
            e = cellfun(@(x) x.exp, results(1,:));
            legends = {e.description};
            plot_results(plots, results, legends, out_folder)
        end
                                                    
    end
    
    if NUMBER_OF_SINGULARITIES.run
        %number_of_singularities(1:10:100, 30)
        number_of_singularities(8:10:100, 30)
    end
    
function number_of_singularities(sing_x, n_repeat)
    global LOG;
    global EPS;
    
    exp_name = 'number_of_singularities';
    exp_description = 'Show a plot of the energy against the number of starting singularities.';
    opt_block_inv_gpu = struct('block_inv_gpu', false, 'CUDAKernel_reduce_cols', false, 'TCODS', true);
    
    [out_folder, log] = get_experiment_folder(exp_name);
    log_bak = LOG;
    LOG = log;
    
    log_and_print(LOG, '%s\r\n', exp_name);
    log_and_print(LOG, '%s\r\n', exp_description);
    
    fp = fullfile('..', 'data', 'bunny.off');
    m = Mesh();
    m.loadTM(fp);
    
    F0 = [1];
    THETA0 = [0];
    V0 = [1, 0, 0];
    N_ITER = 1000;
    DEGREE = 4;
    
    i = 1;
    for n_singularities = sing_x
        log_and_print(LOG, 'n_sing: %d\r\n', n_singularities);
        for j = 1:n_repeat
            log_and_print(LOG, '\tj: %d\r\n', j);
            opt_block_inv_gpu.n_singularities = n_singularities;
            res{i,j} = run_lattice_iter(m, F0, THETA0, DEGREE, N_ITER, opt_block_inv_gpu);
        end
        i = i + 1;
    end
    
    %figure
    %plot(sing_x, cellfun(@(x) x.E, res), '--.', ...
    %    sing_x, cellfun(@(x) sum(abs(x.k)>EPS), res), '--.')
    %legend('Energy', 'Final num of sing', 'Location', 'northwest')
    %xlabel('Number of Singularities')
    %%ylabel('Energy')
    %title('Number of Singularities')
    
    energies = cellfun(@(x) x.E, res);
    singularities = cellfun(@(x) sum(abs(x.k)>EPS), res);
    
    figure
    hold on
    plot(sing_x, mean(energies, 2), '--r', ...
        sing_x, mean(singularities, 2), '--b')
    errorbar(sing_x, mean(energies, 2), std(energies, 0, 2), '.r')
    errorbar(sing_x, mean(singularities, 2), std(singularities, 0, 2), '.b')
    legend('Energy', 'Final num of sing', 'Location', 'northwest')
    hold off
    xlabel('Number of Singularities')
    %ylabel('Energy')
    title('Number of Singularities')
    
    LOG = log_bak;
    
    
function y = get_elapsed_total(x)
    if isfield(x, 'elapsed')
        y = x.elapsed;
    elseif isfield(x, 'elapsed_total')
        y = x.elapsed_total;
    else
        error('Elapsed total field does not exist.')
    end
    
function results = filter_meshes(filepaths, cond_func)
    results = {};
    for i = 1:numel(filepaths)
        fp = filepaths{i};
        mesh = Mesh();
        mesh.loadTM(fp);
        try
            if cond_func(mesh)
                results{end+1} = fp;
                disp(fp)
            end
        catch ex
            fprintf('Caught exception while filtering: %s\r\n', ex.message);
        end
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
    res_iter_gpu_nochol = run_lattice_iter(mesh, ...
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
        res_iterative = run_lattice_iter(mesh, ...
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
        res_iter_gpu = run_lattice_iter(mesh, ...
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
%         res_eigs = run_lattice_iter(mesh, ...
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
        res_iter_gpu_chol = run_lattice_iter(mesh, ...
           F0, ...
           THETA0, ...
           DEGREE, ...
           N_ITER, ...
           opt);
 
        opt.gpu = true;
        opt.chol = false;
        res_iter_gpu_nochol = run_lattice_iter(mesh, ...
            F0, ...
            THETA0, ...
            DEGREE, ...
            N_ITER, ...
            opt);
        
        opt.gpu = false;
        opt.chol = true;
        res_iter_nogpu_chol = run_lattice_iter(mesh, ...
           F0, ...
           THETA0, ...
           DEGREE, ...
           N_ITER, ...
           opt);
        
        opt.gpu = false;
        opt.chol = false;
        res_iter_nogpu_nochol = run_lattice_iter(mesh, ...
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
        res_iter_gpu = run_lattice_iter(mesh, ...
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
            res_iter_gpu_nochol = run_lattice_iter(mesh, ...
                F0, ...
                THETA0, ...
                DEGREE, ...
                N_ITER, ...
                opt);

            opt.gpu = true;
            opt.chol = true;
            res_iter_gpu_chol = run_lattice_iter(mesh, ...
               F0, ...
               THETA0, ...
               DEGREE, ...
               N_ITER, ...
               opt);

            opt.gpu = false;
            opt.chol = true;
            res_iter_nogpu_chol = run_lattice_iter(mesh, ...
               F0, ...
               THETA0, ...
               DEGREE, ...
               N_ITER, ...
               opt);

            opt.gpu = false;
            opt.chol = false;
            res_iter_nogpu_nochol = run_lattice_iter(mesh, ...
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
    %res = nrosy_wrapper(fp, f0, v0, degree);
    %log_and_print(LOG, '%s\r\n', res.result);
    %log_and_print(LOG, 'Status: %d\r\n', res.status);
    res = nrosy_mex(fp, f0, v0, degree);
    log_and_print(LOG, 'E: %g\r\n', res.E);
    
    %elapsed = regexp(res.result, 'Total\s+(\d+(?:\.\d+)?)', 'match');
    %elapsed = regexp(elapsed, '(\d+(?:\.\d+)?)', 'match');
    
    %res.elapsed = str2double(elapsed{:}{:});
    
    res.title = {'MIQ ', ...
        sprintf('$E_{MIQ} = %g$', res.E), ...
        sprintf('Elapsed time: %g', res.elapsed_total)};
    
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
    res_iter_gpu = run_lattice_iter(mesh, ...
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
    
    log_and_print(LOG, 'Testing how high we can go\n');
    %mesh_filenames = {'bunny_res1.off', 'bunny_res2.off', 'bunny_res3.off', 'bunny_res4.off', 'bunny_res5.off', 'bunny_res6.off', 'bunny_res7.off', 'bunny_26k_faces.off', 'bunny_57k_faces.off', 'bunny_64k_faces.off', 'bunny_82k_faces.off', 'bunny_99k_faces.off'};
    %mesh_filenames = {'bunny_26k_faces.off', 'bunny_35k_faces.off', 'bunny_57k_faces.off', 'bunny_64k_faces.off', 'bunny_82k_faces.off', 'bunny_99k_faces.off'};
    mesh_filenames = {'bunny_res1.off', 'bunny_res2.off', 'bunny_res3.off', 'bunny_res4.off', 'bunny_res5.off', 'bunny_res6.off', 'bunny_res7.off', 'bunny_26k_faces.off'};
    r = 1;
    for mname = mesh_filenames
        gd = gpuDevice(1);
        data_folder = find_data_folder();
        fp = fullfile(data_folder, mname{:});
        disp(fp)
    
        mesh = Mesh();
        mesh.loadTM(fp);  
        nV = mesh.nV;
        nF = mesh.nF;
        
        opt.cot = true;
        opt.CUDAKernel = true;
        opt.inv_gpu = true;
        opt.iter_gpu = true;
        
        res_CUDA(r) = run_lattice_iter(mesh, ...
            F0, ...
            THETA0, ...
            DEGREE, ...
            N_ITER, ...
            opt);
        
        opt.inv_gpu = false;
        res_chol(r) = run_lattice_iter(mesh, ...
            F0, ...
            THETA0, ...
            DEGREE, ...
            N_ITER, ...
            opt);
        
        opt.inv_gpu = true;
        opt.CUDAKernel_reduce_cols = true;
        opt.debug_CUDA = false;
        res_CUDA_reduce_cols(r) = run_lattice_iter(mesh, ...
            F0, ...
            THETA0, ...
            DEGREE, ...
            N_ITER, ...
            opt);
        
        n_faces(r) = nF;
        
        r = r + 1;
    end
    
    fig1 = figure();
    
    plot(n_faces, arrayfun(@(x) x.elapsed_total, res_CUDA), '--.', ...
        n_faces, arrayfun(@(x) x.elapsed_total, res_chol), '--.', ...
        n_faces, arrayfun(@(x) x.elapsed_total, res_CUDA_reduce_cols), '--.')
    legend('inv gpu, inner gpu (total)', ...
        'inv chol, inner gpu (total)', ...
        'inv gpu, red cols (total', ...
        'Location', 'northwest')
    title('GPU vs Chol (total)')

    fig2 = figure();
    plot(n_faces, arrayfun(@(x) x.elapsed_setup, res_CUDA), '--.', ...
        n_faces, arrayfun(@(x) x.elapsed_setup, res_chol), '--.', ...
        n_faces, arrayfun(@(x) x.elapsed_setup, res_CUDA_reduce_cols), '--.')
    legend('inv gpu, inner gpu (setup)', ...
        'inv chol, inner gpu (setup)', ...
        'inv gpu, red cols (setup)', ...
        'Location', 'northwest')
    title('GPU vs Chol (setup - dlp, b, k, k0)')
    
    fig3 = figure();
    plot(n_faces, arrayfun(@(x) x.elapsed_inv, res_CUDA), '--.', ...
        n_faces, arrayfun(@(x) x.elapsed_inv, res_chol), '--.', ...
        n_faces, arrayfun(@(x) x.elapsed_inv, res_CUDA_reduce_cols), '--.')
    legend('inv gpu, inner gpu (inv)', ...
        'inv chol, inner gpu (inv)', ...
        'inv gpu, red cols (inv)', ...
        'Location', 'northwest')
    title('GPU vs Chol (inv)')
    
    fig4 = figure();
    plot(n_faces, arrayfun(@(x) x.elapsed_iter, res_CUDA), '--.', ...
        n_faces, arrayfun(@(x) x.elapsed_iter, res_chol), '--.', ...
        n_faces, arrayfun(@(x) x.elapsed_iter, res_CUDA_reduce_cols), '--.')
    legend('inv gpu, inner gpu (loop)', ...
        'inv chol, inner gpu (loop)', ...
        'inv gpu, red cols (loop)', ...
        'Location', 'northwest')
    title('GPU vs Chol (loop)')
    
    fig5 = figure();
    plot(n_faces, arrayfun(@(x) x.elapsed_setup_lap, res_CUDA), '--.', ...
        n_faces, arrayfun(@(x) x.elapsed_setup_lap, res_chol), '--.', ...
        n_faces, arrayfun(@(x) x.elapsed_setup_lap, res_CUDA_reduce_cols), '--.')
    legend('inv gpu, inner gpu (laplacian)', ...
        'inv chol, inner gpu (laplacian)', ...
        'inv gpu, red cols (laplacian)', ...
        'Location', 'northwest')
    title('GPU vs Chol (Laplacian)')
    
    mysave(fig1, 'gpu_vs_chol_total');
    mysave(fig2, 'gpu_vs_chol_setup');
    mysave(fig3, 'gpu_vs_chol_inv');
    mysave(fig4, 'gpu_vs_chol_loop');
    mysave(fig5, 'gpu_vs_chol_lap');
    
    % Plot m and energy
    mesh = Mesh();
    mesh.loadTM('../data/bunny_57k_faces.off');  
    nV = mesh.nV;
    nF = mesh.nF;
        
    opt.cot = true;
    opt.CUDAKernel_reduce_cols = true;
    opt.inv_gpu = false;
    opt.iter_gpu = true;
    opt.plot = false;
    opt.TCODS = false;

    res_plot_cot_lap = run_lattice_iter(mesh, ...
        F0, ...
        THETA0, ...
        DEGREE, ...
        N_ITER, ...
        opt);
    title('bunny.off run hist (cot lap)')
    mysave(gcf, 'bunny_run_hist_cot_lap');
    
    opt.graph = true;
    res_plot_graph_lap = run_lattice_iter(mesh, ...
        F0, ...
        THETA0, ...
        DEGREE, ...
        N_ITER, ...
        opt);
    title('bunny.off run hist (graph lap)')
    mysave(gcf, 'bunny_run_hist_graph_lap');  
%     mesh = Mesh();
%     mesh.loadTM(fp);
%     
%     res_CUDA = run_lattice_iter(mesh, ...
%         F0, ...
%         THETA0, ...
%         DEGREE, ...
%         N_ITER, ...
%         opt);
%     
%     V0 = mesh.V(mesh.F(F0, 2), :) - mesh.V(mesh.F(F0, 1), :);   
%     res_MIQ = run_MIQ(fp, F0, V0, DEGREE);
%     elapsed_MIQ = res_MIQ.elapsed;

function results = merge_results(old_results, new_results, experiments, out_folder)
    assert(size(old_results, 1) == size(new_results, 1))
    assert(size(old_results, 2) == size(new_results, 2))
    assert(size(old_results, 2) == numel(experiments))
    
    n_experiments = numel(experiments);
    
    results = old_results;
    for c = 1:n_experiments
        exp = experiments{c};
        if exp.rerun
            results(:, c) = new_results(:, c);
        end
    end
    
    save(fullfile(out_folder, 'results.mat'), 'results');

function [results, out_folder] = run_experiments(filepaths, experiments, exp_name, description)
    % function [results, out_folder] = run_experiments(filepaths, experiments, exp_name, description)
    % Run the given experiments.
    %
    % Input:
    %   filepaths - list of meshes
    %   experiments - list of experiments. Each experiment is essentially a
    %   function that's to be run with each given mesh. Each entry
    %   in experiments is a structure with the following fields:
    %       func
    %       args
    %       description (optimal)
    %   exp_name - experiment name
    %   description - short description of the experiment (will be printed into LOG)
    %
    % Output:
    %   results - The final results are returned in cellarray as 
    %       results{mesh, exp}
    
    global LOG;
    
    [out_folder, log] = get_experiment_folder(exp_name);
    LOG_BAK = LOG;
    LOG = log;
    
    n_experiments = numel(experiments);
    n_files = numel(filepaths);
    n_experiments_total = n_files * n_experiments;
    
    results = cell(n_files, n_experiments);
    
    ME = [];
    try
        log_and_print(log, '%s\r\n', description);
        gd = gpuDevice();
        % On each mesh
        for r = 1:n_files
            fp = filepaths{r};
            log_and_print(log, 'Mesh: %s\r\n', fp);
            % Run all of the experiments
            for c = 1:n_experiments
                exp = experiments{c};
                if isfield(exp, 'rerun') && ~exp.rerun
                    log_and_print(log, '\tSkipping experiment: %s\r\n', exp.description);
                    results{r, c}.exp = exp;
                    continue
                end
                log_and_print(log, '\tExperiment: %s\r\n', exp.description);
                
                try
                    results{r, c} = feval(exp.func, fp, exp.args{:});
                    results{r, c}.failed = false;
                catch ME
                    log_and_print(log, '\tExperiment failed. \r\n');
                    log_and_print(log, '\tStack: %s\r\n', getReport(ME, 'extended'));
                    results{r, c}.exception = ME;
                    results{r, c}.failed = true;
                end
                % Also save the experiment for future reference
                results{r, c}.exp = exp; 
                reset(gd); wait(gd);
            end
            log_and_print(LOG, '%d // %d\r\n', r*numel(experiments), n_experiments_total);
        end
        
        % Save results
        log_and_print(log, 'Saving results to %s/results.mat...\r\n', out_folder);
        save(fullfile(out_folder, 'results.mat'), 'results');
        log_and_print(log, 'Done.\r\n');
        
    % A hack since matlab does not have a finally clause
    catch ME
    end
    if ~isempty(ME)
        log_and_print(log, 'Some of the experiments failed to run properly. See the log for details: %s/log.txt\r\n', out_folder);
        log_and_print(log, 'Last ME: %s\r\n', getReport(ME, 'extended'));
    end
    % Close open resources
    fclose(log);
    LOG = LOG_BAK;
    
function cpu_inv_VS_MIQ()
    global LOG;
    global DEGREE;
    global OUT_FOLDER;
    
    F0 = [1];
    THETA0 = [0];
    N_ITER = 1000;
    ex = [];
    
    log_and_print(LOG, 'Testing how high we can go\n');
    %mesh_filenames = {'bunny_res1.off', 'bunny_res2.off', 'bunny_res3.off'};
    %mesh_filenames = {'bunny_res1.off', 'bunny_res2.off', 'bunny_res3.off', 'bunny_res4.off', 'bunny_res5.off', 'bunny_res6.off', 'bunny_res7.off', 'bunny_26k_faces.off', 'bunny_57k_faces.off', 'bunny_64k_faces.off', 'bunny_82k_faces.off', 'bunny_99k_faces.off'};
    %mesh_filenames = {'bunny_26k_faces.off', 'bunny_35k_faces.off', 'bunny_57k_faces.off', 'bunny_64k_faces.off', 'bunny_82k_faces.off', 'bunny_99k_faces.off'};
    %mesh_filenames = {'bunny_res1.off', 'bunny_res2.off', 'bunny_res3.off', 'bunny_res4.off', 'bunny_res5.off', 'bunny_res6.off', 'bunny_res7.off', 'bunny_26k_faces.off'};
    mesh_filenames = {'bunny_57k_faces.off'};
    r = 1;
    for mname = mesh_filenames
        gd = gpuDevice(1); wait(gd);
        data_folder = find_data_folder();
        fp = fullfile(data_folder, mname{:});
        disp(fp)
    
        mesh = Mesh();
        mesh.loadTM(fp);  
        nV = mesh.nV;
        nF = mesh.nF;
        n_faces(r) = nF;
        
        opt.cot = true;
        opt.CUDAKernel = true;
        opt.inv_gpu = false;
        opt.iter_gpu = true;
        opt.TCODS = true;
        try
            res_ter = run_lattice_iter(mesh, ...
                F0, ...
                THETA0, ...
                DEGREE, ...
                N_ITER, ...
                opt);
            
            V0 = mesh.V(mesh.F(F0, 2), :) - mesh.V(mesh.F(F0, 1), :);   
            res_MIQ = run_MIQ(fp, F0, V0, DEGREE);
        catch ex
            log_and_print(LOG, 'Failed on %s with %d faces and %d vertices. Will rethrow exception after plotting.\n', fp, nF, nV);
            break
        end
        
        meshes{r} = mesh;
        results_MIQ{r} = res_MIQ;
        results_iter{r} = res_ter;
        
        r = r + 1;  
    end
    
    fig1 = figure();
    plot(n_faces, cellfun(@(x) x.elapsed, results_MIQ), '--.', ...
        n_faces, cellfun(@(x) x.elapsed_total, results_iter), '--.')
    legend('MIQ', ...
        'Iterative (inv chol, loop gpu)', ...
        'Location', 'northwest')
    ylabel('Time (s)');
    xlabel('Number of faces');
    title('MIQ vs invChol+loopGPU timing')

    fig2 = figure();
    plot(n_faces, cellfun(@(x) x.E, results_MIQ), '--.', ...
        n_faces, cellfun(@(x) x.E, results_iter), '--.')
    legend('MIQ', ...
        'Iterative (inv chol, loop gpu)', ...
        'Location', 'northwest')
    ylabel('E');
    xlabel('Number of faces');
    title('MIQ vs invChol+loopGPU energy')

    fig3 = figure();
    plot(n_faces, cellfun(@(x) size(x.S, 1), results_MIQ), '--.', ...
        n_faces, cellfun(@(x) size(x.S, 1), results_iter), '--.')
    legend('MIQ', ...
        'Iterative (gpu, nochol)')
    ylabel('Number of singularities');
    xlabel('Number of faces');
    title('MIQ vs invChol+loopGPU Number of singularities')
    
    out_filename = 'MIQ_vs_invChol_loopGPU';
    mysave(fig1, sprintf('%s_timing', out_filename));
    mysave(fig2, sprintf('%s_energy', out_filename));
    mysave(fig3, sprintf('%s_n_singularities', out_filename));
    
    MeshVis.wfigs([out_filename, '_cross_fields'], ...
                    {meshes{:}, meshes{:}}, ...
                    'OutFolder', OUT_FOLDER, ...
                    'Nrosy', {results_MIQ{:}, results_iter{:}}, ...
                    'ConstrainedFaces', F0, ...
                    'ConstraintVectors', V0, ...
                    'Rows', 2);
                
    if ~isempty(ex)
        rethrow(ex)
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
    

        
           
function str = to_string(x)
    str = evalc(['disp(x)']);
        
    
