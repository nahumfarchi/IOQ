FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = true;
SAVE = false;
LOAD = false;
OUT_FOLDER = 'results';
mkdir(OUT_FOLDER);
FONT_SIZE=12;

fp = '../../../data/rounded_cube_keenan.off';
%fp = '../../../data/ashish_nob/knot1.off';
%fp = '../../../data/torus_fat_r2.off';
%fp = '../../../data/multires/multires/gearbox_multires/mesh_3.off';
%fp = '../../../data/ashish_nob/casting_refined.off';
%fp = '../../../data/ashish_nob/fandisk.off';
[~, meshname, ~] = fileparts(fp);
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;

% run IOQ
rng(SEED)
[alpha, beta, elapsed_ioq] = IOQ_highgenus_gpu(...
    V, F, ...
    'UseGPU', USE_GPU, ...
    'Iterations', 2000);

k = [alpha; beta];
res_ioq = TCODS(m, ...
    'k', k, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);
%[E_edges, Emiq] = per_edge_energy(res_tc);

% run MIQ
[res_miq, elapsed_miq] = ...
    nrosy_mex(fp, FACE0, GVEC, DEGREE);

%% plot
%title_ioq = {'IOQ', ...
%    sprintf('$E_M{theta, p) = %.3f$', res_ioq.miq_energy), ...
%    sprintf('ns = %d', res_ioq.n_vert_sing), ...
%    sprintf('T = %.3f', elapsed_ioq)};
%title_miq = {'MIQ', ...
%    sprintf('$E_M(\\theta, p) = %.3f$', res_miq.miq_energy), ...
%    sprintf('ns = %d', res_miq.n_vert_sing), ...
%    sprintf('T = %.3f', elapsed_miq)};
title_ioq = {'IOQ', ...
    ['$E_M(\theta,p)$=', num2str(res_ioq.miq_energy), ...
    ', ns = ', num2str(res_ioq.n_vert_sing), ...
    ', T = ', num2str(elapsed_ioq)]};
title_miq = {'MIQ', ...
    ['$E_M(\theta,p)$=', num2str(res_miq.miq_energy), ...
    ', ns = ', num2str(res_miq.n_vert_sing), ...
    ', T = ', num2str(elapsed_miq)]};

switch meshname
    case 'torus_fat_r2'
        az = -178.1719;
        el = 0.3135;
    case 'rounded_cube_keenan'
        az = -47.2598;
        el = 26.3027;
    case 'knot1'
        az = 177.4561;
        el = 17.7793;
    case 'casting_refined'
        az = -179.2490;
        el = 0.8477;
    otherwise
        az = 0;
        el = 0;
end
opts = {'EdgeAlpha', 0.2, 'FaceAlpha', 1, 'PlotField', false, 'MarkerSize', 10, 'View', [az, el]};

fh = figure; res_ioq.draw(opts{:});
%view([az, el])
plt = Plot(fh, true);
plt.Title = title_ioq;
plt.FontSize = FONT_SIZE;
%maxfig;

if SAVE
    filename = fullfile('results', [meshname, '_ioq']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure; res_miq.draw(opts{:});
view([az, el])
plt = Plot(fh, true);
plt.Title = title_miq;
plt.FontSize = FONT_SIZE;
%maxfig;
if SAVE
    filename = fullfile('results', [meshname, '_miq']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    %print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end




