
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
USE_GPU = true;
OUT_FOLDER = fullfile('experiments', 'random_init_beta_p', 'results');

%data_folder = '../data/genus1_small';
%filepaths = get_filepaths(data_folder, '.off');
%n_files = numel(filepaths);
fp = '../data/genus1_small/elephant_r.off';
m = Mesh(fp);
V = m.V; F = m.F; ng2 = 2*m.genus; nv = m.nV;
[A, K, d0, d1, H] = tcods_gsystem(V, F);
alpha_G = K(1:end-ng2);
beta_G = K(end-ng2+1:end);
%E = zeros(n_files, N_METHODS);

[~, meshname, ~] = fileparts(fp);
if ~exist(OUT_FOLDER, 'dir')
    mkdir(OUT_FOLDER);
end
out_fp = fullfile(OUT_FOLDER, meshname);


%% try gridsearch2 with random beta_p option 2
UB = 2;
LB = -2;
N_REPEAT = 30;
BETA_P_INIT = zeros(ng2, N_REPEAT);
BETA_P = zeros(ng2, N_REPEAT);
ALPHA_P = zeros(nv, N_REPEAT);
OPT2_N = 10;
%E = zeros(N_REPEAT, OPT2_N);
E_opt2_hists = {};
SEEDS = 1000*(1:N_REPEAT)';
Lp = [];

for i = 1:N_REPEAT
    rng(SEEDS(i));
    bp = round(rand(ng2, 1)*(UB-LB) + LB);
    rng(SEED);
    %rng(SEED);
    if USE_GPU
        [alpha_P,beta_P,Lp,out] = IOQ_highgenus_gpu(...
            V, F, ...
            'Laplacian', LAP_TYPE,...
            'Plot', PLOT,...
            'Iterations', 1000, ...
            'highg_method', 'option2', ...
            'n_alternating', OPT2_N, ...
            'beta_P', bp, ...
            'LaplacianPInv', Lp);
    else
        [alpha_P,beta_P,Lp,E_hist,m_hist, E_outer_hist] = IOQ_highgenus(...
            V, F, ...
            'Laplacian', LAP_TYPE,...
            'Plot', PLOT,...
            'Iterations', 1000, ...
            'highg_method', 'option2', ...
            'n_alternating', OPT2_N, ...
            'beta_P', bp, ...
            'LaplacianPInv', Lp);
    end
    %k = [alpha_P; beta_P];  
    %m = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    %fprintf('E = %g\n', m.miq_energy);
    
    BETA_P_INIT(:, i) = bp;
    BETA_P(:, i) = beta_P;
    ALPHA_P(:, i) = alpha_P;
    %E(i, :) = E_outer_hist;
    E_opt2_hists{end+1} = out;
end

%Eround = zeros(N_REPEAT, OPT2_N);
Eround_opt2_hists = {};
fprintf('round\n')
for i = 1:N_REPEAT
    rng(SEEDS(i));
    if USE_GPU
        [alpha_P,beta_P,Lp,out] = IOQ_highgenus_gpu(...
            V, F, ...
            'Laplacian', LAP_TYPE,...
            'Plot',PLOT,...
            'Iterations',1000, ...
            'highg_method', 'option2', ...
            'n_alternating', OPT2_N, ...
            'beta_P', 'round', ...
            'LaplacianPInv', Lp);
    else
        [alpha_P,beta_P,Lp,E_hist,m_hist,E_outer_hist] = IOQ_highgenus(...
            V, F, ...
            'Laplacian', LAP_TYPE,...
            'Plot', PLOT,...
            'Iterations', 1000, ...
            'highg_method', 'option2', ...
            'n_alternating', OPT2_N, ...
            'beta_P', 'round', ...
            'LaplacianPInv', Lp);
    end
    %Eround(i,:) = E_outer_hist;
    Eround_opt2_hists{end+1} = out;
end

%Ezero = zeros(N_REPEAT, OPT2_N);
Ezero_opt2_hists = {};
fprintf('zero\n')
for i = 1:N_REPEAT
    rng(SEEDS(i));
    if USE_GPU
        [alpha_P,beta_P,Lp,out] = IOQ_highgenus_gpu(...
            V, F, ...
            'Laplacian', LAP_TYPE,...
            'Plot', PLOT,...
            'Iterations', 1000, ...
            'highg_method', 'option2', ...
            'n_alternating', OPT2_N, ...
            'beta_P', [], ...
            'LaplacianPInv', Lp);
    else
        [alpha_P,beta_P,Lp,E_hist,m_hist,E_outer_hist] = IOQ_highgenus(...
            V, F, ...
            'Laplacian', LAP_TYPE,...
            'Plot', PLOT,...
            'Iterations', 1000, ...
            'highg_method', 'option2', ...
            'n_alternating', OPT2_N, ...
            'beta_P', [], ...
            'LaplacianPInv', Lp);
    end
    %Ezero(i,:) = E_outer_hist;
    Ezero_opt2_hists{end+1} = out;
end

close all;

%% Plot
% figure; subplot(121); hold on
% X = E(1:N_REPEAT, :);
% xx = 1:size(X, 2);
% for i = 1:size(X, 2)
%     col = X(:, i);
%     idx = ~isnan(col);
%     yy(i) = mean(col(idx));
%     ee(i) = std(col(idx));
% end
% idx = ~isnan(yy);
% errorbar(xx(idx), yy(idx), ee(idx));
% 
% X = Eround(1:N_REPEAT, :);
% xx = 1:size(X, 2);
% for i = 1:size(X, 2)
%     col = X(:, i);
%     idx = ~isnan(col);
%     yy(i) = mean(col(idx));
%     ee(i) = std(col(idx));
% end
% idx = ~isnan(xx);
% errorbar(xx(idx), yy(idx), ee(idx));
% 
% X = Ezero(1:N_REPEAT, :);
% xx = 1:size(X, 2);
% for i = 1:size(X, 2)
%     col = X(:, i);
%     idx = ~isnan(col);
%     yy(i) = mean(col(idx));
%     ee(i) = std(col(idx));
% end
% idx = ~isnan(xx);
% errorbar(xx(idx), yy(idx), ee(idx));
% hold off
% 
% legend('beta_0=rand([-2,2])', 'beta_0=round(...)', 'beta_0=0')
% title('Energy with random initial beta (option 2)')
% xlabel('Iteration (outer)')
% ylabel('Energy')

%% try gridsearch2 with random beta_p option3_optimized
UB = 2;
LB = -2;
N_REPEAT = 30;
BETA_P_INIT = zeros(ng2, N_REPEAT);
BETA_P = zeros(ng2, N_REPEAT);
ALPHA_P = zeros(nv, N_REPEAT);
OPT2_N = 10;
%E = zeros(N_REPEAT, OPT2_N);
E_opt3_hists = {};
SEEDS = 1000*(1:N_REPEAT)';
Lp = [];

for i = 1:N_REPEAT
    rng(SEEDS(i));
    bp = round(rand(ng2, 1)*(UB-LB) + LB);
    rng(SEED);
    if USE_GPU
        [alpha_P,beta_P,Lp, out] = IOQ_highgenus_gpu(...
            V, F, ...
            'Laplacian', LAP_TYPE,...
            'Plot', PLOT,...
            'Iterations', 1000, ...
            'highg_method', 'option3_optimized', ...
            'n_alternating', OPT2_N, ...
            'beta_P', bp, ...
            'LaplacianPInv', Lp);
    else
        [alpha_P,beta_P,Lp,E_hist,m_hist, E_outer_hist] = IOQ_highgenus(...
            V, F, ...
            'Laplacian', LAP_TYPE,...
            'Plot', PLOT,...
            'Iterations', 1000, ...
            'highg_method', 'option3_optimized', ...
            'n_alternating', OPT2_N, ...
            'beta_P', bp, ...
            'LaplacianPInv', Lp);
    end
    %k = [alpha_P; beta_P];  
    %m = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    %fprintf('E = %g\n', m.miq_energy);
    
    BETA_P_INIT(:, i) = bp;
    BETA_P(:, i) = beta_P;
    ALPHA_P(:, i) = alpha_P;
    E_opt3_hists{end+1} = out;
end

%Eround = zeros(N_REPEAT, OPT2_N);
Eround_opt3_hists = {};
fprintf('round\n')
for i = 1:N_REPEAT
    rng(SEEDS(i));
    if USE_GPU
        [alpha_P,beta_P,Lp,out] = IOQ_highgenus_gpu(...
            V, F, ...
            'Laplacian', LAP_TYPE,...
            'Plot',PLOT,...
            'Iterations',1000, ...
            'highg_method', 'option3_optimized', ...
            'n_alternating', OPT2_N, ...
            'beta_P', 'round', ...
            'LaplacianPInv', Lp);
    else
        [alpha_P,beta_P,Lp,E_hist,m_hist,E_outer_hist] = IOQ_highgenus(...
            V, F, ...
            'Laplacian', LAP_TYPE,...
            'Plot', PLOT,...
            'Iterations', 1000, ...
            'highg_method', 'option3_optimized', ...
            'n_alternating', OPT2_N, ...
            'beta_P', 'round', ...
            'LaplacianPInv', Lp);
    end
    %Eround(i,:) = E_outer_hist;
    Eround_opt3_hists{end+1} = out;
end

%Ezero = zeros(N_REPEAT, OPT2_N);
Ezero_opt3_hists = {};
fprintf('zero\n')
for i = 1:N_REPEAT
    rng(SEEDS(i));
    if USE_GPU
        [alpha_P,beta_P,Lp,out] = IOQ_highgenus_gpu(...
            V, F, ...
            'Laplacian', LAP_TYPE,...
            'Plot', PLOT,...
            'Iterations', 1000, ...
            'highg_method', 'option3_optimized', ...
            'n_alternating', OPT2_N, ...
            'beta_P', [], ...
            'LaplacianPInv', Lp);
    else
        [alpha_P,beta_P,Lp,E_hist,m_hist,E_outer_hist] = IOQ_highgenus(...
            V, F, ...
            'Laplacian', LAP_TYPE,...
            'Plot', PLOT,...
            'Iterations', 1000, ...
            'highg_method', 'option3_optimized', ...
            'n_alternating', OPT2_N, ...
            'beta_P', [], ...
            'LaplacianPInv', Lp);
    end
    %Ezero(i,:) = E_outer_hist;
    Ezero_opt3_hists{end+1} = out;
end

close all;

%% Plot option 2
figure(); 
subplot(221);
y1 = -15; y2 = 70;
ylim([y1, y2])
plot_data(E_opt2_hists);
title('Option 2, mean of 30 trials with random beta_0 in [-2,2]')
xlabel('Iteration')
ylabel('Energy')
xticks(0:2:max(xticks))

subplot(222);
ylim([y1, y2])
plot_data(Eround_opt2_hists);
title('Option 2, mean of 30 trials with beta_0=round(...), alpha_0=rand')
xlabel('Iteration')
ylabel('Energy')
xticks(0:2:max(xticks))

subplot(223)
ylim([y1, y2])
plot_data(Ezero_opt2_hists)
title('Option 2, mean of 30 trials with beta_0=0, alpha_0=rand')
xlabel('Iteration')
ylabel('Energy')
xticks(0:2:max(xticks))

subplot(224)
plot_bars({E_opt2_hists, Eround_opt2_hists, Ezero_opt2_hists}, ...
    {'Emiq', 'E1', 'E2', 'E', 'm'}, ...
    {'Opt 2 rand beta', 'Opt 2 beta=round', 'opt 2 beta=0'});
title('Final energy')
yl = ylim;
ylim([0, yl(2)])

%% plot option 3
figure(); 
subplot(221);
y1 = -15; y2 = 30;
ylim([y1, y2])
plot_data(E_opt3_hists);
title('Option 3, mean of 30 trials with random beta_0 in [-2,2]')
xlabel('Iteration')
ylabel('Energy')
xticks(0:2:max(xticks))

subplot(222);
ylim([y1, y2])
plot_data(Eround_opt3_hists);
title('Option 3, mean of 30 trials with beta_0=round(...), alpha_0=rand')
xlabel('Iteration')
ylabel('Energy')
xticks(0:2:max(xticks))

subplot(223)
ylim([y1, y2])
plot_data(Ezero_opt3_hists)
title('Option 3, mean of 30 trials with beta_0=0, alpha_0=rand')
xlabel('Iteration')
ylabel('Energy')
xticks(0:2:max(xticks))

subplot(224)
plot_bars({E_opt3_hists, Eround_opt3_hists, Ezero_opt3_hists}, ...
    {'Emiq', 'E1', 'E2', 'E', 'm'}, ...
    {'Opt 3 rand beta', 'Opt 3 beta=round', 'opt 3 beta=0'});
title('Final energy')
yl = ylim;
ylim([0, yl(2)])

%% Plot option 2 vs option 3
figure
plot_bars({E_opt2_hists, Eround_opt2_hists, Ezero_opt2_hists, ...
           E_opt3_hists, Eround_opt3_hists, Ezero_opt3_hists}, ...
    {'Emiq', 'E1', 'E2', 'E', 'm'}, ...
    {'Opt 2 rand beta', 'Opt 2 beta=round', 'opt 2 beta=0', ...
     'Opt 3 rand beta', 'Opt 3 beta=round', 'opt 3 beta=0'});
yl = ylim;
ylim([0, yl(2)])
title('Final energy')


%% Save results
save(out_fp, ...
    'E_opt2_hists', 'Eround_opt2_hists', 'Ezero_opt2_hists', ...
    'E_opt3_hists', 'Eround_opt3_hists', 'Ezero_opt3_hists')

% %% Plot
% subplot(122); hold on
% X = E(1:N_REPEAT, :);
% xx = 1:size(X, 2);
% for i = 1:size(X, 2)
%     col = X(:, i);
%     idx = ~isnan(col);
%     yy(i) = mean(col(idx));
%     ee(i) = std(col(idx));
% end
% idx = ~isnan(yy);
% errorbar(xx(idx), yy(idx), ee(idx));
% 
% X = Eround(1:N_REPEAT, :);
% xx = 1:size(X, 2);
% for i = 1:size(X, 2)
%     col = X(:, i);
%     idx = ~isnan(col);
%     yy(i) = mean(col(idx));
%     ee(i) = std(col(idx));
% end
% idx = ~isnan(xx);
% errorbar(xx(idx), yy(idx), ee(idx));
% 
% X = Ezero(1:N_REPEAT, :);
% xx = 1:size(X, 2);
% for i = 1:size(X, 2)
%     col = X(:, i);
%     idx = ~isnan(col);
%     yy(i) = mean(col(idx));
%     ee(i) = std(col(idx));
% end
% idx = ~isnan(xx);
% errorbar(xx(idx), yy(idx), ee(idx));
% hold off
% 
% legend('beta_0=rand([-2,2])', 'beta_0=round(...)', 'beta_0=0')
% title('Energy with random initial beta (option 3 optimized)')
% xlabel('Iteration (outer)')
% ylabel('Energy')
% 
% %% Create latex tables
% col_titles = {};
% for i = 1:OPT2_N
%     col_titles{end+1} = num2str(i);
% end
% row_titles = arrayfun(@(x) 'rand $\beta_p$', ones(N_REPEAT,1), 'UniformOutput', false);
% row_titles = {row_titles{:}, 'round(...)', 'zero'};
% 
% clear input;
% input.tablePlacement = 'H';
% input.tableColumnAlignment = 'X';
% input.min = true;
% input.tableCaption = 'Random initial $\beta_p$';
% 
% latex = bold_latex_table(row_titles, col_titles, E, input);
% latex = {['% ', EXP_DESC], '\afterpage{', '\clearpage', '\thispagestyle{empty}', '\begin{landscape}', EXP_DESC, ...
%     latex{:}, ...
%     '\end{landscape}', '\clearpage', '}'};
% fprintf('\n\n\n\n');
% disp(char(latex))