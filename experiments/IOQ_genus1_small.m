%%
FACE0 = 1;
THETA0 = nan;
GVEC = [1, 0, 0];
DEGREE = 4;
SEED = 112;
%TITLES = {'IOQ\_cpu\_genus0', 'Option 1a', 'run\_lattice\_iter (old version)', 'IOQ\_cpu (miri)', 'MIQ', 'Option 2 (n=1)', 'Option 2 (n=5)', 'Option 2 (n=10)'};
TITLES = {'$\beta_p=0$', 'Opt 1a', 'old', 'IOQ\_cpu (miri)', 'MIQ', 'Opt 2 (n=1)', 'Opt 2 (n=2)', 'Opt 2 (n=3)', 'Optz 2 (n=1)', 'Opt 3', 'Opt 1b'};
N_METHODS = numel(TITLES);
PLOT = false;
LAP_TYPE = 'conn';
LOG = -1;
INIT_BETA_P = 'round';
EXP_DESC = [replace(mfilename, '_', '\_'), ' ', datestr(now), '. Initial $beta_p=round(...)$'];

data_folder = '../data/genus1_small';
filepaths = get_filepaths(data_folder, '.off');
n_files = numel(filepaths);
E = zeros(n_files, N_METHODS);
SINGULARITIES = zeros(n_files, N_METHODS);
TIMINGS = zeros(n_files, N_METHODS);
for i = 1:n_files
    fp = filepaths{i};
    m = Mesh(fp);
    V = m.V; F = m.F; ng2 = 2*m.genus; nv = m.nV;
    [A, K, d0, d1, H] = tcods_gsystem(V, F);
    alpha_G = K(1:end-ng2);
    beta_G = K(end-ng2+1:end);
    print_header(fp);
    
    % Genus 0
    rng(SEED);
    tic
    [alpha_P1,Lp1,E_hist1,m_hist1] = IOQ_genus0_cpu(...
        V, F, ...
        'Laplacian', LAP_TYPE, ...
        'Plot', PLOT, ...
        'Iterations',1000);
    T1 = toc;
    k1 = [alpha_P1; zeros(ng2, 1)];
    %S1v = find(alpha_P1); S1v = [S1v, alpha_P1(S1v)];
    %S1g = [(1:ng2)', zeros(ng2, 1)];

    % Option 1a
    rng(SEED);
    tic
    [alpha_P2,beta_P2,Lp2,E_hist2,m_hist2] = ...
        IOQ_highgenus(V, F, ...
            'Laplacian', LAP_TYPE, ...
            'Plot',PLOT, ...
            'Iterations', 1000, ...
            'highg_method', 'option1a', ...
            'beta_P', INIT_BETA_P);
    T2 = toc;
    k2 = [alpha_P2; beta_P2];
    %S2v = find(alpha_P2); S2v = [S2v, alpha_P2(S2v)];
    %S2g = find(beta_P2); S2g = [S2g, beta_P2(S2g)];
    check_norm('alpha_P1','alpha_P2', 'Log', LOG);

    % Old version
    rng(SEED);
    opt.cot = false;
    tic
    res3=run_lattice_iter(m,FACE0,THETA0,DEGREE,1000,opt);
    T3 = toc;
    k3 = [res3.k; zeros(ng2, 1)];

    % Miri's code
    rng(SEED);
    tic
    [alpha_P4,Lp4,E_hist4] = IOQ_cpu(...
        V, F, ...
        'Laplacian', LAP_TYPE, ...
        'Plot', PLOT);
    T4 = toc;
    a4 = Lp4 * (alpha_G - (pi/2)*alpha_P4);
    beta_P4 = round((beta_G - H'*d0*a4) / (pi/2));
    k4 = [alpha_P4; beta_P4];
    
    % MIQ
    %[theta5, p5, k_miq5, R5, local_frames5, frame_diffs5, elapsed5, pFixed5] = NRoSy_mex_bin(fp, FACE0, GVEC, DEGREE);
    tic
    [theta5, p5, kmiq5, R5, local_frames5, frame_diffs5, elapsed5, pFixed5] = ...
        NRoSy_mex_bin(fp, FACE0, GVEC, DEGREE);
    T5 = toc;
    E5 = E_MIQ(m, theta5, frame_diffs5, p5, DEGREE);
    %alpha_P5 = (4 - (-4*kmiq5)) * (pi / 2);
    
    % Change of vars from MIQ to TCODS
    alpha_P5 = round((DEGREE/(2*pi))*(alpha_G-d0'*frame_diffs5) - d0'*p5);
    assert(check_norm('find(abs(alpha_P5)>1e-10)', 'find(kmiq5)', 'Log', LOG));
    beta_P5 = round(4*((K(end-ng2+1:end) + H'*frame_diffs5) / (2*pi) + 0.25*H'*p5));
    k5 = [alpha_P5; beta_P5];
    %k(abs(k) < 1e-5) = 0;
    x5 = d1'*theta5 + frame_diffs5 + (2*pi/DEGREE)*p5;
    x5 = -x5;
    assert(check_norm('norm(x5)^2', 'E5', 'Log', LOG));
    %S_miq = find(k_miq5(1:mesh.nV)); kk = -k_miq5(S_miq);
    %mesh.vert_sing = [S_miq,kk];
    %figure()
    %mesh.draw()
    
    % Option 2 (n=1)
    rng(SEED);
    tic
    [alpha_P6,beta_P6,Lp6,E_hist6,m_hist6] = ...
        IOQ_highgenus(V, F, ...
            'Laplacian', LAP_TYPE, ...
            'Plot', PLOT, ...
            'Iterations', 1000, ...
            'highg_method', 'option2', ...
            'n_alternating', 1, ...
            'beta_P', INIT_BETA_P);
    T6 = toc;
    k6 = [alpha_P6; beta_P6];
    
    % Option 2 (n=5)
    rng(SEED);
    tic
    [alpha_P7,beta_P7,Lp7,E_hist7,m_hist7] = ...
        IOQ_highgenus(V, F, ...
            'Laplacian', LAP_TYPE, ...
            'Plot',PLOT, ...
            'Iterations', 1000, ...
            'highg_method', 'option2', ...
            'n_alternating', 2, ...
            'beta_P', INIT_BETA_P);
    T7 = toc;
    k7 = [alpha_P7; beta_P7];
    
    % Option 2 (n=10)
    rng(SEED);
    tic
    [alpha_P8,beta_P8,Lp8,E_hist8,m_hist8] = ...
        IOQ_highgenus(V, F, ...
            'Laplacian', LAP_TYPE, ...
            'Plot', PLOT, ...
            'Iterations', 1000, ...
            'highg_method', 'option2', ...
            'n_alternating', 3, ...
            'beta_P', INIT_BETA_P);
    T8 = toc;
    k8 = [alpha_P8; beta_P8];
    
    % Option 2 optimized (n=1)
    rng(SEED);
    tic
    [alpha_P9,beta_P9,Lp9,E_hist9,m_hist9] = ...
        IOQ_highgenus(V, F, ...
            'Laplacian', LAP_TYPE, ...
            'Plot', PLOT, ...
            'Iterations', 1000, ...
            'highg_method', 'option2_optimized', ...
            'n_alternating', 1, ...
            'beta_P', INIT_BETA_P);
    T9 = toc;
    k9 = [alpha_P9; beta_P9];
    
    % Option 3 optimized
    rng(SEED);
    tic
    [alpha_P10,beta_P10,Lp10,E_hist10,m_hist10] = ...
        IOQ_highgenus(V, F, ...
            'Laplacian', LAP_TYPE, ...
            'Plot', PLOT, ...
            'Iterations', 1000, ...
            'highg_method', 'option3_optimized', ...
            'beta_P', INIT_BETA_P);
    T10 = toc;
    k10 = [alpha_P10; beta_P10];
    
    % cvp (opt 1b)
    rng(SEED);
    tic
    [alpha_P11,beta_P11,Lp11,E_hist11,m_hist11] = ...
        IOQ_highgenus(V, F, ...
            'Laplacian', LAP_TYPE, ...
            'Plot', PLOT, ...
            'Iterations', 1000, ...
            'highg_method', 'option1b', ...
            'beta_P', INIT_BETA_P);
    T11 = toc;
    k11 = [alpha_P11; beta_P11];
    
    m1 = TCODS(m, 'k', k1, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    m2 = TCODS(m, 'k', k2, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    m3 = TCODS(m, 'k', k3, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    m4 = TCODS(m, 'k', k4, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    m5 = TCODS(m, 'k', k5, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    %assert(check_norm('m5.miq_energy', 'E5'));
    m6 = TCODS(m, 'k', k6, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    m7 = TCODS(m, 'k', k7, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    m8 = TCODS(m, 'k', k8, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    m9 = TCODS(m, 'k', k9, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    m10 = TCODS(m, 'k', k10, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    m11 = TCODS(m, 'k', k11, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
    fprintf('E1 : %g\nE2 : %g\nE3 : %g\nE4 : %g\nE5 : %g\nE6 : %g\nE7 : %g\nE8 : %g\nE9 : %g\nE10 : %g\nE11 : %g\n\n', m1.miq_energy, m2.miq_energy, m3.miq_energy, m4.miq_energy, E5, m6.miq_energy, m7.miq_energy, m8.miq_energy, m9.miq_energy, m10.miq_energy, m11.miq_energy);
    E(i, 1) = m1.miq_energy;
    E(i, 2) = m2.miq_energy;
    E(i, 3) = m3.miq_energy;
    E(i, 4) = m4.miq_energy;
    E(i, 5) = E5;
    E(i, 6) = m6.miq_energy;
    E(i, 7) = m7.miq_energy;
    E(i, 8) = m8.miq_energy;
    E(i, 9) = m9.miq_energy;
    E(i, 10) = m10.miq_energy;
    E(i, 11) = m11.miq_energy;
    SINGULARITIES(i, 1) = m1.n_singularities;
    SINGULARITIES(i, 2) = m2.n_singularities;
    SINGULARITIES(i, 3) = m3.n_singularities;
    SINGULARITIES(i, 4) = m4.n_singularities;
    SINGULARITIES(i, 5) = nnz(alpha_P5) + nnz(beta_P5);
    SINGULARITIES(i, 6) = m6.n_singularities;
    SINGULARITIES(i, 7) = m7.n_singularities;
    SINGULARITIES(i, 8) = m8.n_singularities;
    SINGULARITIES(i, 9) = m9.n_singularities;
    SINGULARITIES(i, 10) = m10.n_singularities;
    SINGULARITIES(i, 11) = m10.n_singularities;
    TIMINGS(i, 1) = T1;
    TIMINGS(i, 2) = T2;
    TIMINGS(i, 3) = T3;
    TIMINGS(i, 4) = T4;
    TIMINGS(i, 5) = T5;
    TIMINGS(i, 6) = T6;
    TIMINGS(i, 7) = T7;
    TIMINGS(i, 8) = T8;
    TIMINGS(i, 9) = T9;
    TIMINGS(i, 10) = T10;
    TIMINGS(i, 11) = T11;
end

% Create latex tables

meshnames = cell(n_files, 1);
for i = 1:n_files
    fp = filepaths{i};
    [~, mname, ~] = fileparts(fp);
    mname = replace(mname, '_', '\_');
    meshnames{i} = mname;
end

clear input;
input.tablePlacement = 'H';
input.tableColumnAlignment = 'X';
input.min = true;

input.tableCaption = 'MIQ Energy';
[latex_E, T_E] = bold_latex_table(meshnames, TITLES, E, input);
input.tableCaption = 'Timing';
latex_time = bold_latex_table(meshnames, TITLES, TIMINGS, input);
input.tableCaption = 'Singularities';
latex_sing = bold_latex_table(meshnames, TITLES, SINGULARITIES, input);

latex = {['% ', EXP_DESC], '\afterpage{', '\clearpage', '\thispagestyle{empty}', '\begin{landscape}', EXP_DESC, ...
    latex_E{:}, latex_time{:}, latex_sing{:}, ...
    '\end{landscape}', '\clearpage', '}'};
fprintf('\n\n\n\n');
disp(char(latex))


% T = cell(n_files+1, N_METHODS+1);
% T(1, 2:end) = TITLES;
% T(2:end, 1) = meshnames;
% T(2:end, 2:end) = arrayfun(@(x) num2str(x), E, 'UniformOutput', false);
% [~, bold] = min(E, [], 2);
% for i = 1:n_files
%     j = bold(i);
%     T{i+1, j+1} = ['$\textbf{', T{i+1,j+1}, '}$'];
% end
% 
% if isempty(INIT_BETA_P)
%     caption_E = 'MIQ Energy (conn lap). Initial $\beta_p=0$';
%     caption_T = 'Timing (conn lap). Initial $\beta_p=0$';
%     caption_S = 'Singularities (conn lap). Initial $\beta_p=0$';
% elseif ischar(INIT_BETA_P) && strcmpi(INIT_BETA_P, 'round')
%     caption_E = 'MIQ Energy (conn lap). Initial $\beta_p=round(...)$';
%     caption_T = 'Timing (conn lap). Initial $\beta_p=round(...)$';
%     caption_S = 'Singularities (conn lap). Initial $\beta_p=round(...)$';
% else
%     errro('?')
% end
% 
% clear input;
% input.data = table(T);
% input.tableCaption = caption_E;
% input.tablePlacement = 'H';
% input.tableColumnAlignment = 'X';
% latex_E = latexTable(input);
% fprintf('\n\n\n');
% 
% T = cell(n_files+1, N_METHODS+1);
% T(1, 2:end) = TITLES;
% T(2:end, 1) = meshnames;
% T(2:end, 2:end) = arrayfun(@(x) num2str(x), TIMINGS, 'UniformOutput', false);
% [~, bold] = min(TIMINGS, [], 2);
% for i = 1:n_files
%     j = bold(i);
%     T{i+1, j+1} = ['$\textbf{', T{i+1,j+1}, '}$'];
% end
% 
% clear input;
% input.data = table(T);
% input.tableCaption = caption_T;
% input.tablePlacement = 'H';
% input.tableColumnAlignment = 'X';
% latex_time = latexTable(input);
% fprintf('\n\n\n');
% 
% T = cell(n_files+1, N_METHODS+1);
% T(1, 2:end) = TITLES;
% T(2:end, 1) = meshnames;
% T(2:end, 2:end) = arrayfun(@(x) num2str(x), SINGULARITIES, 'UniformOutput', false);
% [~, bold] = min(SINGULARITIES, [], 2);
% for i = 1:n_files
%     j = bold(i);
%     T{i+1, j+1} = ['$\textbf{', T{i+1,j+1}, '}$'];
% end
% 
% clear input;
% input.data = table(T);
% input.tableCaption = caption_S;
% input.tablePlacement = 'H';
% input.tableColumnAlignment = 'X';
% latex_sing = latexTable(input);
% fprintf('\n\n\n');


%% try gridsearch2 with random beta_p
% UB = 10;
% LB = -10;
% N_ITER = 10;
% BETA_P_INIT = zeros(ng2, N_ITER);
% BETA_P = zeros(ng2, N_ITER);
% ALPHA_P = zeros(nv, N_ITER);
% E = zeros(N_ITER, 1);
% OPT2_N = 20;
% SEEDS = 1000*(1:N_ITER)';
% for i = 1:N_ITER
%     rng(SEEDS(i));
%     bp = round(rand(ng2, 1)*(UB-LB) + LB);
%     rng(SEED);
%     [alpha_P,beta_P,Lp,E_hist,m_hist] = IOQ_highgenus(...
%         V, F, ...
%         'Laplacian', LAP_TYPE,...
%         'Plot',PLOT,...
%         'Iterations',1000, ...
%         'highg_method', 'option2', ...
%         'n_alternating', OPT2_N, ...
%         'beta_P', bp);
%     k = [alpha_P; beta_P];  
%     m = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE);
%     fprintf('E = %g\n', m.miq_energy);
%     
%     BETA_P_INIT(:, i) = bp;
%     BETA_P(:, i) = beta_P;
%     ALPHA_P(:, i) = alpha_P;
%     E(i) = m.miq_energy;
% end
% 
% fprintf('round\n')
% rng(SEED);
% [alpha_P,beta_P,Lp,E_hist,m_hist] = IOQ_highgenus(...
%     V, F, ...
%     'Laplacian', LAP_TYPE,...
%     'Plot',PLOT,...
%     'Iterations',1000, ...
%     'highg_method', 'option2', ...
%     'n_alternating', OPT2_N, ...
%     'beta_P', 'round');
% k = [alpha_P; beta_P];  
% m = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE);
% fprintf('E = %g\n', m.miq_energy);
% E_round = m.miq_energy;
% 
% fprintf('zero\n')
% rng(SEED);
% [alpha_P,beta_P,Lp,E_hist,m_hist] = IOQ_highgenus(...
%     V, F, ...
%     'Laplacian', LAP_TYPE,...
%     'Plot',PLOT,...
%     'Iterations',1000, ...
%     'highg_method', 'option2', ...
%     'n_alternating', OPT2_N, ...
%     'beta_P', []);
% k = [alpha_P; beta_P];  
% m = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE);
% fprintf('E = %g\n', m.miq_energy);
% E_zero = m.miq_energy;
% 
% E
% E_round
% E_zero
% 
% BETA_P_INIT
% BETA_P