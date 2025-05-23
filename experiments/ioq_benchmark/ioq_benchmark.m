ioq_benchmark_setup;
filepaths = get_filepaths(data_folder, data_ext);

set(groot, 'defaultAxesTickLabelInterpreter','none'); 
set(groot, 'defaultLegendInterpreter','none');
set(0,'defaulttextinterpreter','none')

log_and_print(LOG, create_header(['experiment_name : ', experiment_name]));
log_and_print('Results are being saved into %s\r\n', base_folder);
log_and_print('log file: %s\r\n', fullfile(base_folder, 'log.txt'));

BLOCK_SIZE = 20000;
experiments = {
    struct(...
 		'func', @miq_wrapper, ...
 		'args', {{F0, G_CONSTRAINT_VEC, DEGREE}}, ...
 		'description', 'Mixed integer quadrangulation', ...
 		'name', 'MIQ', ...
 		'save_ffield', true, ...
 		'save_nrosy', true, ...
 		'rerun', true, ...
 		'legend', 'MIQ'), ...
    struct(...
 		'func', @go_wrapper, ...
 		'args', {{DEGREE}}, ...
 		'description', 'Globally optimal direction fields', ...
 		'name', 'GO', ...
 		'save_ffield', true, ...
 		'save_nrosy', true, ...
 		'rerun', true, ...
 		'legend', 'GO'), ...
    struct(...
        'func', @ioq_wrapper, ...
        'args', {{'Laplacian', 'conn', ...
                  'highg_method', 'option1a', ...
                  'beta_P', 'round', ...
                  'BlockSize', BLOCK_SIZE}}, ...
        'description', 'IOQ with connectivity laplacian', ...
        'name', 'IOQ_conn_1a', ...
        'save_ffield', true, ...
        'save_nrosy', true, ...
        'rerun', true, ...
        'legend', 'IOQ (conn, opt1a)'), ...
    struct(...
        'func', @ioq_wrapper, ...
        'args', {{'Laplacian', 'conn', ...
                  'highg_method', 'option1a', ...
                  'beta_P', 'round', ...
                  'BlockSize', BLOCK_SIZE, ...
                  'InvMethod', 'ApproxResistance', ...
                  'JLEps', 0.3, ...
                  'Iterations', 2000}}, ...
        'description', 'JL IOQ with eps=0.3, conn lap', ...
        'name', 'JL_IOQ_eps0.3', ...
        'save_ffield', true, ...
        'save_nrosy', true, ...
        'rerun', true, ...
        'legend', 'JL IOQ (eps=0.3)'), ...
    struct(...
        'func', @ioq_wrapper, ...
        'args', {{'Laplacian', 'conn', ...
                  'highg_method', 'option1a', ...
                  'beta_P', 'round', ...
                  'BlockSize', BLOCK_SIZE, ...
                  'InvMethod', 'ApproxResistance', ...
                  'JLEps', 0.5, ...
                  'Iterations', 2000}}, ...
        'description', 'JL IOQ with eps=0.5, conn lap', ...
        'name', 'JL_IOQ_eps0.5', ...
        'save_ffield', true, ...
        'save_nrosy', true, ...
        'rerun', true, ...
        'legend', 'JL IOQ (eps=0.5)'), ...
    struct(...
        'func', @ioq_wrapper, ...
        'args', {{'Laplacian', 'conn', ...
                  'highg_method', 'option1a', ...
                  'beta_P', 'round', ...
                  'BlockSize', BLOCK_SIZE, ...
                  'InvMethod', 'ApproxResistance', ...
                  'JLEps', 0.8, ...
                  'Iterations', 2000}}, ...
        'description', 'JL IOQ with eps=0.8, conn lap', ...
        'name', 'JL_IOQ_eps0.8', ...
        'save_ffield', true, ...
        'save_nrosy', true, ...
        'rerun', true, ...
        'legend', 'JL IOQ (eps=0.8)'), ...
    };
    
% experiments = {
% 	struct(...
% 		'func', @miq_wrapper, ...
% 		'args', {{F0, G_CONSTRAINT_VEC, DEGREE}}, ...
% 		'description', 'MIQ', ...
% 		'name', 'MIQ', ...
% 		'save_ffield', true, ...
% 		'save_nrosy', true, ...
% 		'rerun', false, ...
% 		'legend', 'MIQ'), ...
%     struct(...
%         'func', @ioq_wrapper, ...
%         'args', {{'Laplacian', 'conn', ...
%                   'highg_method', 'genus0', ...
%                   'beta_P', 'round', ...
%                   'BlockSize', BLOCK_SIZE}}, ...
%         'description', 'IOQ with connectivity laplacian, genus0', ...
%         'name', 'IOQ_conn_genus0', ...
%         'save_ffield', true, ...
%         'save_nrosy', true, ...
%         'rerun', false, ...
%         'legend', 'IOQ (conn, genus0)'), ...
%     struct(...
%         'func', @ioq_wrapper, ...
%         'args', {{'Laplacian', 'cot', ...
%                   'highg_method', 'genus0', ...
%                   'beta_P', 'round', ...
%                   'BlockSize', BLOCK_SIZE}}, ...
%         'description', 'IOQ with cot laplacian, genus0', ...
%         'name', 'IOQ_cot_genus0', ...
%         'save_ffield', true, ...
%         'save_nrosy', true, ...
%         'rerun', false, ...
%         'legend', 'IOQ (cot, genus0)'), ...
%     struct(...
%         'func', @ioq_wrapper, ...
%         'args', {{'Laplacian', 'conn', ...
%                   'highg_method', 'option1a', ...
%                   'beta_P', 'round', ...
%                   'BlockSize', BLOCK_SIZE}}, ...
%         'description', 'IOQ with connectivity laplacian, opt1a', ...
%         'name', 'IOQ_conn_opt1a', ...
%         'save_ffield', true, ...
%         'save_nrosy', true, ...
%         'rerun', false, ...
%         'legend', 'IOQ (conn, opt1a)'), ...
%     struct(...
%         'func', @ioq_wrapper, ...
%         'args', {{'Laplacian', 'cot', ...
%                   'highg_method', 'option1a', ...
%                   'beta_P', 'round', ...
%                   'BlockSize', BLOCK_SIZE}}, ...
%         'description', 'IOQ with cot laplacian, opt1a', ...
%         'name', 'IOQ_cot_opt1a', ...
%         'save_ffield', true, ...
%         'save_nrosy', true, ...
%         'rerun', false, ...
%         'legend', 'IOQ (cot, opt1a)'), ...
%     struct(...
%         'func', @ioq_wrapper, ...
%         'args', {{'Laplacian', 'conn', ...
%                   'highg_method', 'option2_optimized', ...
%                   'beta_P', 'round', ...
%                   'BlockSize', BLOCK_SIZE}}, ...
%         'description', 'IOQ with connectivity laplacian, opt2z', ...
%         'name', 'IOQ_conn_opt2z', ...
%         'save_ffield', true, ...
%         'save_nrosy', true, ...
%         'rerun', false, ...
%         'legend', 'IOQ (conn, opt2z)'), ...
%     struct(...
%         'func', @ioq_wrapper, ...
%         'args', {{'Laplacian', 'cot', ...
%                   'highg_method', 'option2_optimized', ...
%                   'beta_P', 'round', ...
%                   'BlockSize', BLOCK_SIZE}}, ...
%         'description', 'IOQ with cot laplacian, opt2z', ...
%         'name', 'IOQ_cot_opt2z', ...
%         'save_ffield', true, ...
%         'save_nrosy', true, ...
%         'rerun', false, ...
%         'legend', 'IOQ (cot, opt2z)'), ...
%     struct(...
%         'func', @ioq_wrapper, ...
%         'args', {{'Laplacian', 'conn', ...
%                   'highg_method', 'option1a', ...
%                   'beta_P', 'round', ...
%                   'BlockSize', BLOCK_SIZE, ...
%                   'InvMethod', 'ApproxResistance', ...
%                   'JLEps', 0.5, ...
%                   'Iterations', 2000}}, ...
%         'description', 'JL IOQ with eps=0.5, conn lap', ...
%         'name', 'JL_IOQ_eps0.5', ...
%         'save_ffield', true, ...
%         'save_nrosy', true, ...
%         'rerun', true, ...
%         'legend', 'JL IOQ (eps=0.5)'), ...
%     struct(...
%         'func', @ioq_wrapper, ...
%         'args', {{'Laplacian', 'conn', ...
%                   'highg_method', 'option1a', ...
%                   'beta_P', 'round', ...
%                   'BlockSize', BLOCK_SIZE, ...
%                   'InvMethod', 'ApproxResistance', ...
%                   'JLEps', 0.8, ...
%                   'Iterations', 2000}}, ...
%         'description', 'JL IOQ with eps=0.8, conn lap', ...
%         'name', 'JL_IOQ_eps0.8', ...
%         'save_ffield', true, ...
%         'save_nrosy', true, ...
%         'rerun', false, ...
%         'legend', 'JL IOQ (eps=0.8)'), ...
%     struct(...
%         'func', @ioq_wrapper, ...
%         'args', {{'Laplacian', 'conn', ...
%                   'highg_method', 'option1a', ...
%                   'beta_P', 'round', ...
%                   'BlockSize', BLOCK_SIZE, ...
%                   'InvMethod', 'ApproxResistance', ...
%                   'JLEps', 1, ...
%                   'Iterations', 2000}}, ...
%         'description', 'JL IOQ with eps=1, conn lap', ...
%         'name', 'JL_IOQ_eps1', ...
%         'save_ffield', true, ...
%         'save_nrosy', true, ...
%         'rerun', false, ...
%         'legend', 'JL IOQ (eps=1)'), ...
% };

n_experiments = numel(experiments);
n_files = numel(filepaths);
SEED = 112;

%stats_fp = fullfile(...
%    '..', 'results', 'experiments', experiment_name, 'stats.mat');
stats_fp = fullfile(base_folder, 'stats.mat');
if exist(stats_fp, 'file')
    stats = load(stats_fp);
    stats = stats.stats;
else
    stats = cell(n_files+1, n_experiments+1);
end

%% Create and save ffields
counter = 1;
progressbar('Files', 'Experiments')
failures = {};
for r = 1:n_files
	fp = filepaths{r};
    %fp = '../../../data/ashish_nob_small/block.off';
	log_and_print(LOG, 'Mesh: %s\r\n', fp);

	[~, name, ~] = fileparts(fp);
	stats{r+1, 1} = name;

	energies = [];

	for c = 1:n_experiments
        %print_header(sprintf('%d / %d', counter, n_files*n_experiments));
        log_and_print(LOG, create_header(sprintf('%d / %d', counter, n_files*n_experiments)));
		exp = experiments{c};
		stats{1, c+1} = exp.name;

		log_and_print(LOG, '\tExperiment: %s\r\n', exp.description);
        if isfield(exp, 'rerun') && ~exp.rerun
            log_and_print(LOG, '\tSkipping...\r\n');
            counter = counter + 1;
            continue
        end
        
        [~, fname, ~] = fileparts(fp);
        out_fp_ffield = fullfile(out_folder_ffields, ...
            [fname, '_', exp.name, ffield_ext]);
        %[~, fname, ~] = fileparts(fp);
	    out_fp_nrosy = fullfile(out_folder_nrosy, ...
            [fname, '_', exp.name, nrosy_ext]);
        
        %if isfield(stats{r+1, c+1}, 'failed') && ~stats{r+1, c+1}.failed && ...
        %        exist(out_fp_ffield, 'file') && ...
        %        exist(out_fp_nrosy, 'file')
        %    log_and_print(LOG, '\tSkipping...\r\n');
        %    counter = counter + 1;
        %    continue
        %end
        
        try
            rng(SEED);
            res = feval(exp.func, fp, exp.args{:});                

	        if isfield(exp, 'save_ffield') && exp.save_ffield
	            log_and_print(LOG, 'Saving %s...\r\n', out_fp_ffield);
	            res.m.saveFField(out_fp_ffield);	            
	        end
	        if isfield(exp, 'save_nrosy') && exp.save_nrosy
	            log_and_print(LOG, 'Saving %s...\r\n', out_fp_nrosy);
	            res.m.saveNRosy(out_fp_nrosy, true);
	        end
	        
	        st.nV = res.m.nV;
	        st.nF = res.m.nF;
	        st.nE = res.m.nE;
            st.genus = res.m.genus;
	        st.miq_energy = res.m.miq_energy;
            st.miq_edge = per_edge_energy(res.m);
	        st.n_singularities = res.m.n_singularities;
	        st.elapsed_total = res.elapsed_total;
	        st.exp = exp; 
	        st.failed = false;
	        st.alg_used = res.alg_used;
	        %st.symbol = SYMBOLS(res.alg_used);
	        st.fp = fp;
	        stats{r+1, c+1} = st;
	        
	        energies(end+1) = res.m.miq_energy;	
        catch ME
            log_and_print(LOG, '\tExperiment failed. \r\n');
            log_and_print(LOG, '\tStack: %s\r\n', getReport(ME, 'extended'));
            stats{r+1, c+1}.failed = true;
            stats{r+1, c+1}.exp = exp;
            failures{end+1} = exp;
        end
        
        counter = counter + 1;
        save(stats_fp, 'stats');
        
        frac2 = c / n_experiments;
        frac1 = ((r-1) + frac2) / n_files;
        progressbar(frac1, frac2)
	end

	disp(['Energies (', name, ')'])
    disp(energies)
end
progressbar(1, 1)

%% Save statistics
log_and_print(LOG, 'Saving statistics to %s...\r\n', stats_fp);
save(stats_fp, 'stats');
log_and_print(LOG, 'Done.\r\n');