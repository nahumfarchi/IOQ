FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = true;
SAVE = false;
LOAD = false;
OUT_FOLDER = 'results';
mkdir(OUT_FOLDER);
FONT_SIZE=12;
COLORS = linspecer(3);
COLORS = {COLORS(1, :), COLORS(3, :), COLORS(2, :)};

fp = '../../../data/bunny.off';
%fp = '../../../data/elephant_r.off';
%fp = '../../../data/basic/cat.off';
[~, meshname, ~] = fileparts(fp);
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;

N_INIT_SING = [0, 50, 100];
out_exact = {};
for i = 1:length(N_INIT_SING)
    ns = N_INIT_SING(i);
    rng(SEED)
    [alpha, beta, elapsed_ioq, ~, out1] = IOQ_highgenus_gpu(...
        V, F, ...
        'UseGPU', USE_GPU, ...
        'Iterations', 2000, ...
        'Histories', true, ...
        'NSingularities', ns);
    out_exact{i} = out1;
end

out_approx = {};
for i = 1:length(N_INIT_SING)
    ns = N_INIT_SING(i);
    rng(SEED)
    [alpha, beta, elapsed_ioq, ~, out1] = IOQ_highgenus_gpu(...
        V, F, ...
        'UseGPU', USE_GPU, ...
        'Iterations', 2000, ...
        'Histories', true, ...
        'NSingularities', ns, ...
        'InvMethod', 'ApproxResistance', ...
        'JLEps', 1);
    out_approx{i} = out1;
end

%%
close all
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex')

%fh = figure;
tt1_e = 1:length(out_exact{1}.E_hist1);
tt2_e = 1:length(out_exact{2}.E_hist1);
tt3_e = 1:length(out_exact{3}.E_hist1);
yy1_e = out_exact{1}.E_hist1;
yy2_e = out_exact{2}.E_hist1;
yy3_e = out_exact{3}.E_hist1;
plt = Plot(tt1_e, yy1_e, tt2_e, yy2_e, tt3_e, yy3_e);
plt.ShowBox = 'off';
%plt = Plot(fh, true);
plt.Colors = COLORS;
plt.Legend = arrayfun(@(x) num2str(x), N_INIT_SING, 'UniformOutput', false);
plt.LegendBox = 'on';
plt.XLabel = 'Iteration';
plt.YLabel = 'Energy';
plt.XLim = [0, length(tt3_e)];

export_fig('convergence.pdf');
export_fig('convergence.png');

tt1_a = 1:length(out_approx{1}.E_hist1);
tt2_a = 1:length(out_approx{2}.E_hist1);
tt3_a = 1:length(out_approx{3}.E_hist1);
yy1_a = out_approx{1}.E_hist1;
yy2_a = out_approx{2}.E_hist1;
yy3_a = out_approx{3}.E_hist1;
plt = Plot(tt1_a, yy1_a, tt2_a, yy2_a, tt3_a, yy3_a);
plt.ShowBox = 'off';
%plt = Plot(fh, true);
plt.Colors = COLORS;
plt.Legend = arrayfun(@(x) num2str(x), N_INIT_SING, 'UniformOutput', false);
plt.LegendBox = 'on';
plt.XLabel = 'Iteration';
plt.YLabel = 'Energy';
plt.XLim = [0, length(tt3_a)];

export_fig('convergence_approx.pdf');
export_fig('convergence_approx.png');

%%
COLORS = linspecer(6);
COLORS = COLORS(end:-1:1, :);
plt = Plot(...
    tt1_e, yy1_e, tt2_e, yy2_e, tt3_e, yy3_e, ...
    tt1_a, yy1_a, tt2_a, yy2_a, tt3_a, yy3_a);
plt.Colors = {COLORS(1,:), COLORS(2,:), COLORS(3,:), COLORS(4,:), COLORS(5,:), COLORS(6,:)};
plt.Legend = arrayfun(@(x) num2str(x), [N_INIT_SING, N_INIT_SING], 'UniformOutput', false);
plt.LegendBox = 'on';
plt.XLabel = 'Iteration';
plt.YLabel = 'Energy';
plt.LineWidth = 2.5*ones(6, 1);
plt.LineStyle = {'-', '-', '-', ':', ':', ':'}
plt.Legend = {'IOQ $|S|=0$', 'IOQ $|S|=50$', 'IOQ $|S|=100$', 'IOQe $|S|=0$', 'IOQe $|S|=50$', 'IOQe $|S|=100$'};
plt.XLim = [0, max(length(tt3_e), length(tt3_a))];
plt.ShowBox = 'off';

export_fig('convergence_both.pdf')
export_fig('convergence_both.png')

legend('off');
plt.XLim = [44.6815   63.6433];
plt.YLim = [3.9461   37.3956];

export_fig('convergence_both_closeup.pdf')
export_fig('convergence_both_closeup.png')

% plt = Plot(tt3_e, yy3_e, tt3_a, yy3_a);
% plt.Colors = {COLORS(1,:), COLORS(2,:), COLORS(3,:), COLORS(4,:), COLORS(5,:), COLORS(6,:)};
% plt.Legend = arrayfun(@(x) num2str(x), [N_INIT_SING, N_INIT_SING], 'UniformOutput', false);
% plt.LegendBox = 'on';
% plt.XLabel = 'Iteration';
% plt.YLabel = 'Energy';
% plt.LineWidth = ones(6, 1);
% plt.LineStyle = {'-', '-', '-', '--', '--', '--'}
% plt.Legend = {'IOQ ns=100', 'IOQe ns=100'};

%%
set(groot, 'defaultAxesTickLabelInterpreter','none'); 
set(groot, 'defaultLegendInterpreter','none');
set(0,'defaulttextinterpreter','none')

