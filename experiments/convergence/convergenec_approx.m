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
[~, meshname, ~] = fileparts(fp);
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;
d0 = get_exterior_derivatives(m);
L = d0'*d0;
Lp = invChol_mex(full(L + 1/nv)) - 1/nv;

N_INIT_SING = [0, 50, 100];
out = {};
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
        'LaplacianPInv', Lp);
    out{i} = out1;
end

%%
close all
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex')

%fh = figure;
tt1 = 1:length(out{1}.E_hist1);
tt2 = 1:length(out{2}.E_hist1);
tt3 = 1:length(out{3}.E_hist1);
plt = Plot(tt1, out{1}.E_hist1, tt2, out{2}.E_hist1, tt3, out{3}.E_hist1);
plt.ShowBox = 'off';
%plt = Plot(fh, true);
plt.Colors = COLORS;
plt.Legend = arrayfun(@(x) num2str(x), N_INIT_SING, 'UniformOutput', false);
plt.LegendBox = 'on';
plt.XLabel = 'Iteration';
plt.YLabel = 'Energy';

export_fig('convergence_approx.pdf');
export_fig('convergence_approx.png');

set(groot, 'defaultAxesTickLabelInterpreter','none'); 
set(groot, 'defaultLegendInterpreter','none');
set(0,'defaulttextinterpreter','none')

