%% ========================================================================
%% Test and time various Lx=b solvers for JL IOQ.
%% ========================================================================

FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0];

data_folder = '../../../data/bunnies_small';
[~, n, ~] = fileparts(data_folder);
out_folder = fullfile('results', n);
if ~exist(out_folder, 'dir')
    mkdir(out_folder);
end
filepaths = get_filepaths(data_folder);
n_files = numel(filepaths);

SAVE = false;
REPS = 1;
SAME_SEED = true;
SEED = 112;
USE_GPU = true;
EPS = 0.5;

E = [];
T = [];
N_FACES = [];
names = {};
rng(SEED);
counter = 1;
for i = n_files:n_files
    fp = filepaths{i};
    print_header(fp);
    m = Mesh(fp); V = m.V; F = m.F; ne = m.nE;
    
    for r = 1:REPS
        fprintf('%d // %g\n', counter, n_files*REPS);
        
        % -----------------------------------------------------------------
        % Approx resistance, 
        % Setup: ichol + backslash
        % Gridsearch: regular
        % -----------------------------------------------------------------
        if SAME_SEED, rng(SEED); end
        [alpha_p, beta_p, elapsed_total] = IOQ_highgenus_gpu(V, F, ...
                    'InvMethod', 'ApproxResistance', ...
                    'highg_method', 'genus0', ...
                    'Iterations', 30, ...
                    'UseGPU', USE_GPU, ...
                    'Mesh', m, ...
                    'bsx', false, ...
                    'JLEps', EPS, ...
                    'Colamd', false);
        k = [alpha_p; beta_p];
        res = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
        E(end+1) = res.miq_energy;
        T(end+1) = elapsed_total;
        names{end+1} = 'ichol+backslash';
        
        % -----------------------------------------------------------------
        % Approx resistance, 
        % Setup: ichol + 1 backslash
        % Gridsearch: regular
        % -----------------------------------------------------------------
        if SAME_SEED, rng(SEED); end
        [alpha_p, beta_p, elapsed_total] = IOQ_highgenus_gpu(V, F, ...
                    'InvMethod', 'ApproxResistance', ...
                    'highg_method', 'genus0', ...
                    'Iterations', 30, ...
                    'UseGPU', USE_GPU, ...
                    'Mesh', m, ...
                    'bsx', false, ...
                    'JLEps', EPS, ...
                    'Colamd', false, ...
                    'ZtildeOneBS', true);
        k = [alpha_p; beta_p];
        res = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
        E(end+1) = res.miq_energy;
        T(end+1) = elapsed_total;
        names{end+1} = 'ichol+1backslash';
                
        % -----------------------------------------------------------------
        % Approx resistance
        % Setup: colamd ichol + backslash
        % Gridsearch: regular
        % -----------------------------------------------------------------
        if SAME_SEED, rng(SEED); end
        [alpha_p, beta_p, elapsed_total] = IOQ_highgenus_gpu(V, F, ...
                    'InvMethod', 'ApproxResistance', ...
                    'highg_method', 'genus0', ...
                    'Iterations', 30, ...
                    'UseGPU', USE_GPU, ...
                    'Mesh', m, ...
                    'bsx', false, ...
                    'JLEps', EPS, ...
                    'Colamd', true);
        k = [alpha_p; beta_p];
        res = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
        E(end+1) = res.miq_energy;
        T(end+1) = elapsed_total;
        names{end+1} = 'coldamd_ichol+backslash';
        
        % -----------------------------------------------------------------
        % Approx resistance
        % Setup: ichol with gpu + backslash
        % Gridsearch: regular
        % -----------------------------------------------------------------
        if SAME_SEED, rng(SEED); end
        [alpha_p, beta_p, elapsed_total] = IOQ_highgenus_gpu(V, F, ...
                    'InvMethod', 'ApproxResistance', ...
                    'highg_method', 'genus0', ...
                    'Iterations', 30, ...
                    'UseGPU', USE_GPU, ...
                    'Mesh', m, ...
                    'bsx', false, ...
                    'JLEps', EPS, ...
                    'Colamd', true, ...
                    'CholGPU', false, ...
                    'ZtildeOneBS', true);
        k = [alpha_p; beta_p];
        res = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
        E(end+1) = res.miq_energy;
        T(end+1) = elapsed_total;
        names{end+1} = 'ichol_gpu+backslash';
    end
    
end

results = table(E', T', 'RowNames', names, 'VariableNames', {'Energy', 'Time'})