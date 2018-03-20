%%
FACE0 = 1;
THETA0 = nan;
GVEC = [1, 0, 0];
DEGREE = 4;
SEED = 112;
%TITLES = {'IOQ\_cpu\_genus0', 'Option 1a', 'run\_lattice\_iter (old version)', 'IOQ\_cpu (miri)', 'MIQ', 'Option 2 (n=1)', 'Option 2 (n=5)', 'Option 2 (n=10)'};
TITLES = {'$\beta_p=0$', 'Option 1a', 'old version', 'IOQ\_cpu (miri)', 'MIQ', 'Option 2 (n=1)', 'Option 2 (n=5)', 'Option 2 (n=10)'};
N_METHODS = numel(TITLES);
PLOT = true;
LAP_TYPE = 'conn';
LOG = -1;
INIT_BETA_P = [];
EXP_DESC = [replace(mfilename, '_', '\_'), ' ', datestr(now), '. Initial random $\beta_p \in [-2,2]$'];

%data_folder = '../data/genus1_small';
%filepaths = get_filepaths(data_folder, '.off');
%n_files = numel(filepaths);
fp = '../data/3holes.off';
m = Mesh(fp);
V = m.V; F = m.F; ng2 = 2*m.genus; nv = m.nV;
[A, K, d0, d1, H] = tcods_gsystem(V, F);
alpha_G = K(1:end-ng2);
beta_G = K(end-ng2+1:end);
%E = zeros(n_files, N_METHODS);


%% try gridsearch2 with random beta_p
UB = 2;
LB = -2;
N_ITER = 30;
BETA_P_INIT = zeros(ng2, N_ITER);
BETA_P = zeros(ng2, N_ITER);
ALPHA_P = zeros(nv, N_ITER);
OPT2_N = 10;
E = zeros(N_ITER, OPT2_N);
SEEDS = 1000*(1:N_ITER)';
rng(SEED);
bp = round(rand(ng2, 1)*(UB-LB) + LB);
for i = 1:N_ITER
    rng(SEEDS(i));
    %rng(SEED);
    [alpha_P,beta_P,Lp,E_hist,m_hist, E_outer_hist] = IOQ_highgenus(...
        V, F, ...
        'Laplacian', LAP_TYPE,...
        'Plot',PLOT,...
        'Iterations',1000, ...
        'highg_method', 'option2', ...
        'n_alternating', OPT2_N, ...
        'beta_P', bp);
    %k = [alpha_P; beta_P];  
    %m = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    %fprintf('E = %g\n', m.miq_energy);
    
    BETA_P_INIT(:, i) = bp;
    BETA_P(:, i) = beta_P;
    ALPHA_P(:, i) = alpha_P;
    E(i, :) = E_outer_hist;
end

Eround = zeros(N_ITER, OPT2_N);
fprintf('round\n')
for i = 1:N_ITER
    rng(SEEDS(i));
    [alpha_P,beta_P,Lp,E_hist,m_hist,E_outer_hist] = IOQ_highgenus(...
        V, F, ...
        'Laplacian', LAP_TYPE,...
        'Plot',PLOT,...
        'Iterations',1000, ...
        'highg_method', 'option2', ...
        'n_alternating', OPT2_N, ...
        'beta_P', 'round');
    Eround(i,:) = E_outer_hist;
end

Ezero = zeros(N_ITER, OPT2_N);
fprintf('zero\n')
for i = 1:N_ITER
    rng(SEEDS(i));
    [alpha_P,beta_P,Lp,E_hist,m_hist,E_outer_hist] = IOQ_highgenus(...
        V, F, ...
        'Laplacian', LAP_TYPE,...
        'Plot',PLOT,...
        'Iterations',1000, ...
        'highg_method', 'option2', ...
        'n_alternating', OPT2_N, ...
        'beta_P', []);
    Ezero(i,:) = E_outer_hist;
end

%% Plot
figure; hold on
X = E(1:N_ITER, :);
xx = 1:size(X, 2);
yy = mean(X, 1);
ee = std(X, 1);
errorbar(xx, yy, ee);

X = Eround(1:N_ITER, :);
xx = 1:size(X, 2);
yy = mean(X, 1);
ee = std(X, 1);
errorbar(xx, yy, ee);

X = Ezero(1:N_ITER, :);
xx = 1:size(X, 2);
yy = mean(X, 1);
ee = std(X, 1);
errorbar(xx, yy, ee);
hold off

legend('beta_0=rand([-2,2])', 'beta_0=round(...)', 'beta_0=0')
title('Energy with random initial beta')




%% Create latex tables
col_titles = {};
for i = 1:OPT2_N
    col_titles{end+1} = num2str(i);
end
row_titles = arrayfun(@(x) 'rand $\beta_p$', ones(N_ITER,1), 'UniformOutput', false);
row_titles = {row_titles{:}, 'round(...)', 'zero'};

clear input;
input.tablePlacement = 'H';
input.tableColumnAlignment = 'X';
input.min = true;
input.tableCaption = 'Random initial $\beta_p$';

latex = bold_latex_table(row_titles, col_titles, E, input);
latex = {['% ', EXP_DESC], '\afterpage{', '\clearpage', '\thispagestyle{empty}', '\begin{landscape}', EXP_DESC, ...
    latex{:}, ...
    '\end{landscape}', '\clearpage', '}'};
fprintf('\n\n\n\n');
disp(char(latex))