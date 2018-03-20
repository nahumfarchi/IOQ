% Create and save ffields from the meshes in data_folder.
% Selects the best possible combinations of cpu/blockgpu/gpu inv and cpu/gpu loop for each mesh.
%
% Results are saved as .ffield and .rosy files.
% See Mesh.saveFField and Mesh.saveNRosy for a description of these formats.
%
% You can load an .ffield file by using:
%   m = Mesh();
%   m.loadFField('bunny.ffield');
%   m.draw();

%%
experiments_setup;
filepaths = get_filepaths(data_folder, data_ext);

% Create algorithm selector (can take a while...)
%log_and_print(LOG, 'Checking ranges for the various inversion methods...\r\n');
%tic
%[alg_selector, ...
%    inv_gpu_range, inv_block_range, inv_cpu_range, iter_gpu_range] = find_IOQ_ranges(filepaths, LB);
%toc
% Note that the ranges are with respect to nV.
%
% On cggc-miri-x299, you can use:
inv_gpu_range = [-inf, 21582];
inv_block_range = [21582, 50002];
inv_cpu_range = [50002, inf];
iter_gpu_range = [-inf, 21582];
alg_selector = create_IOQ_alg_selector(...
        inv_gpu_range, ...
        inv_block_range, ...
        inv_cpu_range, ...
        iter_gpu_range);

%%

print_header(['experiment_name : ', experiment_name])

disp(['Results are being saved into ', base_folder])

opt_conn = struct('TCODS', true, ...
    'gConstraintVec', G_CONSTRAINT_VEC, ...
    'graph', true, ...
    'alg_selector', alg_selector, ...
    'lb', LB);

opt_cot = struct('TCODS', true, ...
    'gConstraintVec', G_CONSTRAINT_VEC, ...
    'graph', false, ...
    'alg_selector', alg_selector, ...
    'lb', LB);

experiments = {
    struct(...
        'func', @run_lattice_iter, ...
        'args', {{F0, THETA0, DEGREE, N_ITER, opt_conn}}, ...
        'description', 'IOQ with connectivity lap', ...
        'name', 'IOQ_conn', ...
        'save_ffield', true, ...
        'save_nrosy', true, ...
        'rerun', true, ...
        'legend', 'IOQ (conn)'), ...
    struct(...
        'func', @run_lattice_iter, ...
        'args', {{F0, THETA0, DEGREE, N_ITER, opt_cot}}, ...
        'description', 'IOQ with cot lap', ...
        'name', 'IOQ_cot', ...
        'save_ffield', true, ...
        'save_nrosy', true, ...
        'rerun', true, ...
        'legend', 'IOQ (cot)'), ...
    struct(...
        'func', @nrosy_mex, ...
        'args', {{F0, G_CONSTRAINT_VEC, DEGREE}}, ...
        'description', 'MIQ', ...
        'name', 'MIQ', ...
        'save_ffield', true, ...
        'save_nrosy', true, ...
        'rerun', false, ...
        'legend', 'MIQ'), ...
};

%% Create and save ffields
n_files = numel(filepaths);
n_experiments = numel(experiments);

rand_seed = rng;

for r = 1:n_files
    fp = filepaths{r};
    log_and_print(LOG, 'Mesh: %s\r\n', fp);
    
    [~, name, ~] = fileparts(fp);
    stats{r+1, 1} = name;
    
    energies = [];
    
    for c = 1:n_experiments
        exp = experiments{c};
        stats{1, c+1} = exp.name;
        
        if isfield(exp, 'rerun') && ~exp.rerun
            log_and_print(LOG, '\tSkipping experiment: %s\r\n', exp.description);
            continue;
        end
        
        log_and_print(LOG, '\tExperiment: %s\r\n', exp.description);
        log_and_print(LOG, '\t%g%%...\r\n', ((r-1)*n_experiments+c-1) / (n_files*n_experiments));
        
        try
            % Reset the seed and run the experiment
            rng(rand_seed);
            res = feval(exp.func, fp, exp.args{:});

            % Save results
            m = Mesh();
            m.loadTM(fp);
            m.set_ffield(...
                res.degree, ...
                res.local_frames, ...
                res.frame_diffs, ...
                res.theta, ...
                res.ffield, ...
                res.S, ...
                res.E)

            if isfield(exp, 'save_ffield') && exp.save_ffield
                [~, fname, ~] = fileparts(fp);
                out_fp_ffield = fullfile(out_folder_ffields, ...
                    [fname, '_', exp.name, ffield_ext]);
                log_and_print(LOG, 'Saving %s...\r\n', out_fp_ffield);
                m.saveFField(out_fp_ffield);
                
            end
            if isfield(exp, 'save_nrosy') && exp.save_nrosy
                [~, fname, ~] = fileparts(fp);
                out_fp_nrosy = fullfile(out_folder_nrosy, ...
                    [fname, '_', exp.name, nrosy_ext]);
                m.saveNRosy(out_fp_nrosy, true);
            end
            
            st.nV = m.nV;
            st.nF = m.nF;
            st.nE = m.nE;
            st.miq_energy = m.miq_energy;
            st.n_singularities = m.n_singularities;
            st.elapsed_total = res.elapsed_total;
            st.exp = exp; 
            st.failed = false;
            st.alg_used = res.alg_used;
            st.symbol = SYMBOLS(res.alg_used);
            st.fp = fp;
            stats{r+1, c+1} = st;
            
            energies(end+1) = m.miq_energy;
        catch ME
            log_and_print(LOG, '\tExperiment failed. \r\n');
            log_and_print(LOG, '\tStack: %s\r\n', getReport(ME, 'extended'));
            stats{r+1, c+1}.failed = true;
            stats{r+1, c+1}.exp = exp;
        end
    
    end
    
    disp(['Energies (', name, ')'])
    disp(energies)
end

%% Save statistics
log_and_print(LOG, 'Saving statistics to %s...\r\n', stats_fp);
save(stats_fp, 'stats');
log_and_print(LOG, 'Done.\r\n');

% %% Plots
% plots = {...
%     struct(...
%         'title', 'Energy', ...
%         'name', 'energy', ...
%         'get_x_field', @(x) x.nF, ...
%         'get_y_field', @(x) x.miq_energy, ...
%         'style', {{'x'}}), ...
%     struct(...
%         'title', 'Elapsed time (total)', ...
%         'name', 'elapsed_total', ...
%         'get_x_field', @(x) x.nF, ...
%         'get_y_field', @(x) x.elapsed_total, ...
%         'style', {{'x'}}), ...
%     struct(...
%         'title', 'Number of Singularities', ...
%         'name', 'n_singularities', ...
%         'get_x_field', @(x) x.nF, ...
%         'get_y_field', @(x) x.n_singularities, ...
%         'style', {{'x'}})
%     };
% 
% e = cellfun(@(x) x.exp, stats(1,:));
% legends = {e.description};
% plot_results(plots, stats, legends, out_folder_plots)
% 
% 
%% Plot IOQ (graph) vs IOQ (cot) vs MIQ
% base_folder = fullfile('..', 'results', 'experiments');
% exp_folder_graph = fullfile(base_folder, 'genus0_connectivity_lap');
% exp_folder_cot = fullfile(base_folder, 'genus0_cot_lap');
% stats_graph = load(fullfile(exp_folder_graph, 'stats.mat'));
% stats_cot = load(fullfile(exp_folder_cot, 'stats.mat'));
% stats_graph = stats_graph.stats;
% stats_cot = stats_cot.stats;
% 
% miq_stats = stats_graph(:, 5);
% IOQ_graph_stats = stats_graph(:, 1:4);
% IOQ_cot_stats = stats_cot(:, 1:4);
% 
% xx = cellfun(@(x) x.nF, miq_stats);
% [xx, inds] = sort(xx);
% 
% miq_stats = miq_stats(inds);
% IOQ_graph_stats = IOQ_graph_stats(inds, :);
% IOQ_cot_stats = IOQ_cot_stats(inds, :);
% 
% % Elapsed time
% % connection laplacian
% figure
% is_valid = @(x) ~isempty(x) && ...
%                 (~isfield(x, 'failed') || ...
%                 ~x.failed);
% valid = cellfun(@(x) is_valid(x), IOQ_graph_stats);
% [rows, cols] = size(IOQ_graph_stats);
% yy_IOQ_graph = inf(rows, cols);
% yy_IOQ_graph(valid) = cellfun(@(x) x.elapsed_total, IOQ_graph_stats(valid));
% [yy_IOQ_graph, I_graph] = min(yy_IOQ_graph, [], 2);
% %yy_IOQ_graph = yy_IOQ_graph(inds);
% %I_graph = I_graph(inds);
% 
% % cot laplacian
% valid = cellfun(@(x) is_valid(x), IOQ_cot_stats);
% [rows, cols] = size(IOQ_cot_stats);
% yy_IOQ_cot = inf(rows, cols);
% yy_IOQ_cot(valid) = cellfun(@(x) x.elapsed_total, IOQ_cot_stats(valid));
% [yy_IOQ_cot, I_cot] = min(yy_IOQ_cot, [], 2);
% %yy_IOQ_cot = yy_IOQ_cot(inds);
% %I_cot = I_cot(inds);
% 
% % miq
% yy_miq = cellfun(@(x) x.elapsed_total, miq_stats);
% %yy_miq = yy_miq(inds);
% 
% plot(xx, yy_IOQ_graph, '--', ...
%     xx, yy_IOQ_cot, '--', ...
%     xx, yy_miq, ['--', experiments{5}.symbol]);
% legends = {'IOQ (connectivity)', 'IOQ (cot)', 'MIQ'};
% 
% symbol_legends = {'cpuinv, cpuloop', ...
%     'cpuinv, gpuloop', ...
%     'blockinv, gpuloop', ...
%     'gpuinv, gpuloop'};
% hold on
% for c = 1:4
%     idx = I_graph == c;
%     if nnz(idx) > 0
%         scatter(xx(idx), yy_IOQ_graph(idx), experiments{c}.symbol);
%         legends{end+1} = symbol_legends{c};
%     end
% end
% hold off
% 
% hold on
% for c = 1:4
%     idx = I_cot == c;
%     if nnz(idx) > 0
%         scatter(xx(idx), yy_IOQ_cot(idx), experiments{c}.symbol);
%     end
% end
% hold off
% 
% legend(legends, 'Location', 'northwest');
% title('Elapsed time')
% 
% % Energy
% % TODO
% figure
% yy_IOQ_graph = inf(rows, cols);
% yy_IOQ_graph(valid) = cellfun(@(x) x.miq_energy, IOQ_graph_stats(valid));
% I = (1 : rows)';
% J = I_graph(:);
% k = sub2ind([rows, cols], I, J);
% yy_IOQ_graph = yy_IOQ_graph(k);
% 
% yy_IOQ_cot = inf(rows, cols);
% yy_IOQ_cot(valid) = cellfun(@(x) x.miq_energy, IOQ_cot_stats(valid));
% I = (1 : rows)';
% J = I_cot(:);
% k = sub2ind([rows, cols], I, J);
% yy_IOQ_cot = yy_IOQ_cot(k);
% 
% yy_miq = cellfun(@(x) x.miq_energy, miq_stats);
% 
% plot(xx, yy_IOQ_graph, '--', ...
%     xx, yy_IOQ_cot, '--', ...
%     xx, yy_miq, ['--', experiments{5}.symbol]);
% legends = {'IOQ (connectivity)', 'IOQ (cot)', 'MIQ'};
% 
% hold on
% for c = 1:4
%     idx = I_graph == c;
%     if nnz(idx) > 0
%         scatter(xx(idx), yy_IOQ_graph(idx), experiments{c}.symbol);
%         legends{end+1} = symbol_legends{c};
%     end
% end
% hold off
% 
% hold on
% for c = 1:4
%     idx = I_cot == c;
%     if nnz(idx) > 0
%         scatter(xx(idx), yy_IOQ_cot(idx), experiments{c}.symbol);
%     end
% end
% hold off
% 
% legend(legends, 'Location', 'northwest');
% title('energy')

% Number of singularities
% TODO
% 
% 
% %%
% plots = {...
%     struct(...
%         'title', 'Elapsed Time', ...
%         'name', 'elapsed_time', ...
%         'get_x_field', @(x) x.nF, ...
%         'get_y_field', @(x) x.miq_energy, ...
%         'style', {{'--'}}, ...
%         'x_col', 5, ...
%         'y_cols', 1:4, ...
%         'symbols', {{'*', 'X', 'O', '^', 'd'}}, ...
%         'symbol_legends', {{'cpuinv, cpuloop', 'cpuinv, gpuloop', 'blockinv, gpuloop', 'gpuinv, gpuloop', 'MIQ'}}, ...
%         'legend', 'IOQ (graph)'), ...
%     struct(...
%         'title', 'Energy', ...
%         'name', 'energy', ...
%         'get_x_field', @(x) x.nF, ...
%         'get_y_field', @(x) x.miq_energy, ...
%         'style', {{'--'}}, ...
%         'x_col', 5, ...
%         'y_cols', 5, ...
%         'symbols', {{'*', 'X', 'O', '^', 'd'}}, ...
%         'symbol_legends', {{'cpuinv, cpuloop', 'cpuinv, gpuloop', 'blockinv, gpuloop', 'gpuinv, gpuloop', 'MIQ'}}, ...
%         'legend', 'MIQ'), ...
%         };
% figure()
% plot_results(plots, stats_graph, {}, []);