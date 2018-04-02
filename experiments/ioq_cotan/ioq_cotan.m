FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = true;
SAVE = false;
LOAD = false;
OUT_FOLDER = 'results';
mkdir(OUT_FOLDER);
FONT_SIZE=14;

%fp = '../../../data/rounded_cube_keenan.off';
%fp = '../../../data/ashish_nob/knot1.off';
%fp = '../../../data/torus_fat_r2.off';
%fp = '../../../data/multires/multires/gearbox_multires/mesh_3.off';
%fp = '../../../data/ashish_nob/casting_refined.off';
%fp = '../../../data/ashish_nob/fandisk.off';
%fp = '../../../data/bunny.off';
%fp = '../../../data/torus_fat_r2.off';
%fp = '../../../data/phands.off';
%fp = '../../../data/round_cuber.off';
%fp = '../../../data/ashish_nob/armadillo.off';
%fp = '../../../data/ashish_nob/robocat_deci.off';
%fp = '../../../data/bunny.off';
%fp = '../../../data/phands_r.off';
%fp = '../../../data/bunny_r.off';
%fp = '../../../data/fandisk_r.off';
fp = '../../../data/bumpy.off';

[~, meshname, ~] = fileparts(fp);
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;

% 
% figure
% subplot(121); m.draw(props{:});
% subplot(122); m_flipped.draw(props{:});

%% run IOQ with connectivity lap
rng(SEED)
[alpha1, beta1, elapsed_ioq1] = IOQ_highgenus_gpu(...
    V, F, ...
    'UseGPU', USE_GPU, ...
    'Iterations', 2000, ...
    'Laplacian', 'conn');

%
k1 = [alpha1; beta1];
res_ioq_conn = TCODS(m, ...
    'k', k1, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);
%[E_edges, Emiq] = per_edge_energy(res_tc);

%% run IOQ with cot lap
rng(SEED)
[alpha2, beta2, elapsed_ioq2] = IOQ_highgenus_gpu(...
    V, F, ...
    'UseGPU', USE_GPU, ...
    'Iterations', 2000, ...
    'Laplacian', 'cot');

%
k2 = [alpha2; beta2];
res_ioq_cot = TCODS(m, ...
    'k', k2, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);
%[E_edges, Emiq] = per_edge_energy(res_tc);

%%
% run MIQ
[res_miq, elapsed_miq] = ...
    nrosy_mex(fp, FACE0, GVEC, DEGREE);

%% Plots
props = {'FaceColor', 'w', 'EdgeColor', 'k'};
title1 = {'IOQ', sprintf('IOQ, $E = %g, |S| = %d$',   res_ioq_conn.miq_energy, res_ioq_conn.n_vert_sing)};
title2 = {'IOQm', sprintf('IOQm, $E = %g, |S| = %d$', res_ioq_cot.miq_energy,  res_ioq_cot.n_vert_sing)};
title3 = {'MIQ', sprintf('MIQ, $E = %g, |S| = %d$',   res_miq.miq_energy,      res_miq.n_vert_sing)};

%% Compare with MIQ
% figure
% subplot(131); res_ioq_conn.draw(props{:}) ; title(title1)
% subplot(132); res_ioq_cot.draw(props{:})  ; title(title2)
% subplot(133); res_miq.draw(props{:})      ; title(title3)

%
figure
subplot(131); res_ioq_conn.draw(props{:})    ; title(title1)
subplot(132); res_ioq_cot.draw(props{:})     ; title(title2)
subplot(133); res_miq.draw(props{:})         ; title(title3)

%set(gcf, 'WindowStyle', 'docked')
%filename = sprintf('%s_%d_flips.png', meshname, n_flips);
%export_fig(filename);
%filename = sprintf('%s_%d_flips.fig', meshname, n_flips);
%savefig(filename)

%%
%props = {'FaceColor', 'w', 'PlotField', true};
%figure
%subplot(121); res_ioq_conn_f.draw(props{:}) ; title(title4)
%subplot(122); res_ioq_cot_f.draw(props{:})  ; title(title5)

% %% plot
% %title_ioq = {'IOQ', ...
% %    sprintf('$E_M{theta, p) = %.3f$', res_ioq.miq_energy), ...
% %    sprintf('ns = %d', res_ioq.n_vert_sing), ...
% %    sprintf('T = %.3f', elapsed_ioq)};
% %title_miq = {'MIQ', ...
% %    sprintf('$E_M(\\theta, p) = %.3f$', res_miq.miq_energy), ...
% %    sprintf('ns = %d', res_miq.n_vert_sing), ...
% %    sprintf('T = %.3f', elapsed_miq)};
% title_ioq = {'IOQ', ...
%     ['$E_M(\theta,p)$=', num2str(res_ioq.miq_energy), ...
%     ', ns = ', num2str(res_ioq.n_vert_sing), ...
%     ', T = ', num2str(elapsed_ioq1)]};
% title_miq = {'MIQ', ...
%     ['$E_M(\theta,p)$=', num2str(res_miq.miq_energy), ...
%     ', ns = ', num2str(res_miq.n_vert_sing), ...
%     ', T = ', num2str(elapsed_miq)]};
% 
% switch meshname
%     case 'torus_fat_r2'
%         az = -178.1719;
%         el = 0.3135;
%     case 'rounded_cube_keenan'
%         az = -47.2598;
%         el = 26.3027;
%     case 'knot1'
%         az = 177.4561;
%         el = 17.7793;
%     case 'casting_refined'
%         az = -179.2490;
%         el = 0.8477;
%     otherwise
%         az = 0;
%         el = 0;
% end
% opts = {'EdgeAlpha', 0.2, 'FaceAlpha', 1, 'PlotField', false, 'MarkerSize', 10, 'View', [az, el]};
% 
% fh = figure; res_ioq.draw(opts{:});
% %view([az, el])
% plt = Plot(fh, true);
% plt.Title = title_ioq;
% plt.FontSize = FONT_SIZE;
% %maxfig;
% 
% if SAVE
%     filename = fullfile('results', [meshname, '_ioq']);
%     export_fig([filename '.pdf'])
%     export_fig([filename '.png'])
%     saveas(gcf, filename, 'fig')
% end
% 
% fh = figure; res_miq.draw(opts{:});
% view([az, el])
% plt = Plot(fh, true);
% plt.Title = title_miq;
% plt.FontSize = FONT_SIZE;
% %maxfig;
% if SAVE
%     filename = fullfile('results', [meshname, '_miq']);
%     export_fig([filename '.pdf'])
%     export_fig([filename '.png'])
%     %print(gcf, filename, '-dpng', '-r300')
%     saveas(gcf, filename, 'fig')
% end




