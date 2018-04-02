FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = true;
SAVE = false;
LOAD = false;
FONT_SIZE = 14;
INV_METHOD = 'GPUInv';

%fp = '../../../data/rounded_cube_keenan.off';
%fp = '../../../data/ashish_nob/knot1.off';
%fp = '../../../data/torus_fat_r2.off';
%fp = '../../../data/multires/multires/gearbox_multires/mesh_3.off';
%fp = '../../../data/ashish_nob/casting_refined.off';
%fp = '../../../data/ashish_nob/fandisk.off';
fp = '../../../data/decimated-max.off';
%fp = '../../../data/torus_fat_r2.off';
%fp = '../../../data/phands.off';
%fp = '../../../data/round_cuber.off';
%fp = '../../../data/ashish_nob/armadillo.off';
%fp = '../../../data/ashish_nob/robocat_deci.off';
%fp = '../../../data/bunny.off';
%fp = '../../../data/phands_r.off';
%fp = '../../../data/bunny_r.off';
%fp = '../../../data/fandisk_r.off';
%fp = '../../../data/bumpy.off';

[~, meshname, ~] = fileparts(fp);
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;

out_folder = fullfile('results', meshname);
mkdir(out_folder);

n_flips = 20;
rng(SEED)
[F_flipped, n_flipped] = flip_edges(m.V, m.F, n_flips, pi/10);
m_flipped = Mesh(m.V, F_flipped);
m_flipped.saveTM(fullfile(out_folder, sprintf('%s_%d.off', meshname, n_flips)));

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
    'Laplacian', 'conn', ...
    'InvMethod', INV_METHOD);
a1 = L \ ((pi/2)*alpha1 - alpha_g);

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
    'Laplacian', 'cot', ...
    'InvMethod', INV_METHOD);

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

%% Run with flipped edges
%% IOQ with conn lap
rng(SEED)
[alpha3, beta3, elapsed_ioq3, ~, out3] = IOQ_highgenus_gpu(...
    m_flipped.V, m_flipped.F, ...
    'UseGPU', USE_GPU, ...
    'Iterations', 2000, ...
    'Laplacian', 'conn', ...
    'InvMethod', INV_METHOD);

%
k3 = [alpha3; beta3];
res_ioq_conn_f = TCODS(m_flipped, ...
    'k', k3, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);
%[E_edges, Emiq] = per_edge_energy(res_tc);

%% run IOQ with cot lap
rng(SEED)
[alpha4, beta4, elapsed_ioq4, Lp, out4] = IOQ_highgenus_gpu(...
    m_flipped.V, m_flipped.F, ...
    'UseGPU', USE_GPU, ...
    'Iterations', 3000, ...
    'Laplacian', 'cot', ...
    'Histories', true, ...
    'InvMethod', INV_METHOD);

%
k4 = [alpha4; beta4];
res_ioq_cot_f = TCODS(m_flipped, ...
    'k', k4, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);

%%
% run MIQ
% [res_miq, elapsed_miq] = ...
%     nrosy_mex(fp, FACE0, GVEC, DEGREE);
% 
% m_flipped.saveTM('tmp.off');
% [res_miq_flipped, elapsed_miq] = ...
%     nrosy_mex('tmp.off', FACE0, GVEC, DEGREE);


%% Plots
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex')

title_ioqe = {'IOQe', sprintf('$E = %g, |S| = %d$',           res_ioq_conn.miq_energy,   res_ioq_conn.n_vert_sing)};
title_ioqe_cot = {'IOQe cot', sprintf('$E = %g, |S| = %d$',   res_ioq_cot.miq_energy,    res_ioq_cot.n_vert_sing)};
title_ioqe_f = {'IOQe', sprintf('$E = %g, |S| = %d$',         res_ioq_conn_f.miq_energy, res_ioq_conn_f.n_vert_sing)};
title_ioqe_cot_f = {'IOQe cot', sprintf('$E = %g, |S| = %d$', res_ioq_cot_f.miq_energy,  res_ioq_cot_f.n_vert_sing)};

% phands camera
switch meshname
    case "phands_r"
        %th.Position=[0.66,0,-0.7]
        cam.pba = [1.307915 1.000000 1.241718 ];
        cam.dar = [1 1 1 ];
        cam.cva = [8.045990 ];
        cam.cuv = [-0.991047 0.130604 -0.027718 ];
        cam.ct = [0.091304 -0.002161 0.003669 ];
        cam.cp = [0.406029 0.785555 -7.537573 ];
    otherwise
        cam = [];
end
%%
CM = cbrewer('div', 'RdBu', 256);
%CM = [];
props = {'FaceAlpha', 1, ...
       'EdgeAlpha', 1, ...
       'Dock', true, ...
       'Colormap', CM, ...
       'Caxis', 'auto', ...
       'Camera', cam, ...
       'EdgeColor', 'none', ...
       'MarkerSize', 40};

% Compare with MIQ
% figure
% subplot(131); res_ioq_conn.draw(props{:}) ; title(title1)
% subplot(132); res_ioq_cot.draw(props{:})  ; title(title2)
% subplot(133); res_miq.draw(props{:})      ; title(title3)

% Compare with edge flips
figure
subplot(221); res_ioq_conn.draw(props{:})    ; title(title_ioqe);
subplot(222); res_ioq_cot.draw(props{:})     ; title(title_ioqe_cot) ;
subplot(223); res_ioq_conn_f.draw(props{:})  ; title(title_ioqe_f) ; 
subplot(224); res_ioq_cot_f.draw(props{:})   ; title(title_ioqe_cot_f) ;

return

%%

FM = 2;
names = {};
SPC = '';
POS = [0.8, 0, -0.7];

figure
filename = fullfile(out_folder, [meshname, '_01.pdf']);
names{end+1} = filename;
res_ioq_conn.draw(props{:})
th = title([title_ioqe, SPC], 'FontSize', FM*FONT_SIZE);
th.Position=POS;
export_fig(filename)
%filename = fullfile(out_folder, [meshname, '_title_01.pdf']);
%create_title(filename, title_ioqe, FONT_SIZE);

figure
filename = fullfile(out_folder, [meshname, '_02.pdf']);
names{end+1} = filename;
res_ioq_conn_f.draw(props{:})
th = title([title_ioqe_f, SPC], 'FontSize', FM*FONT_SIZE);
th.Position=POS;
export_fig(filename)
%filename = fullfile(out_folder, [meshname, '_title_02.pdf']);
%create_title(filename, title_ioqe_f, FONT_SIZE)

figure
filename = fullfile(out_folder, [meshname, '_03.pdf']);
names{end+1} = filename;
res_miq.draw(props{:})
th = title([title_miq, SPC], 'FontSize', FM*FONT_SIZE);
th.Position=POS;
export_fig(filename)
%filename = fullfile(out_folder, [meshname, '_title_03.pdf']);
%create_title(filename, title_miq, FONT_SIZE)

figure
filename = fullfile(out_folder, [meshname, '_04.pdf']);
names{end+1} = filename;
res_miq_flipped.draw(props{:})
th = title(title_miq_f, 'FontSize', FM*FONT_SIZE);
th.Position=POS;
export_fig(filename)
%filename = fullfile(out_folder, [meshname, '_title_04.pdf']);
%create_title(filename, title_miq_f, FONT_SIZE)

out_pan = fullfile(out_folder, 'edge_flips.pdf');
panorama(names, out_pan);
%cmd = ['"C:/Program Files (x86)/IrfanView/i_view32.exe" /panorama=(1'];
%for j=1:length(names)
%    cmd = [cmd ',' names{j}];
%end
%cmd = [cmd ') /convert=' out_pan];
%system(cmd)



return

%%
%set(gcf, 'WindowStyle', 'docked')
filename = sprintf('results/%s_%d_flips.png', meshname, n_flips);
export_fig(filename);
filename = sprintf('results/%s_%d_flips.fig', meshname, n_flips);
savefig(filename)

%% Compare with edge flips (cot)
figure
subplot(231); res_ioq_conn.draw(props{:})    ; title(title1);
subplot(232); res_ioq_cot.draw(props{:})     ; title(title2) ;
subplot(233); res_miq.draw(props{:})         ; title(title3) ; 
subplot(234); res_ioq_conn_f.draw(props{:})  ; title(title4) ; 
subplot(235); res_ioq_cot_f.draw(props{:})   ; title(title5) ;
subplot(236); res_miq_flipped.draw(props{:}) ; title(title6) ; 

% phands back view
switch meshname
    case "phands_r"
        cam.pba = [1.316639 1.006670 1.000000 ];
        cam.dar = [1 1 1 ];
        cam.cva = [6.608610 ];
        cam.cuv = [-0.976437 0.052212 -0.209393 ];
        cam.ct = [0.037156 0.001752 0.001818 ];
        cam.cp = [-1.633013 -2.052472 7.277905 ];
    case "fandisk_r"
        cam.pba = [119.507246 143.288372 119.507246 ];
        cam.dar = [1 1 1 ];
        cam.cva = [6.608610 ];
        cam.cuv = [-0.213468 -0.794841 -0.568031 ];
        cam.ct = [2.413970 15.240700 -1.327014 ];
        cam.cp = [14.112631 -11.913513 32.273199 ];
    otherwise
        cam = [];
end
props = {'FaceAlpha', 1, ...
       'Dock', true, ...
       'Colormap', CM, ...
       'Caxis', 'auto', ...
       'Camera', cam, ...
       'EdgeColor', 'k', ...
       'MarkerSize', 40};

figure
subplot(231); res_ioq_conn.draw(props{:})    ; title(title1);
subplot(232); res_ioq_cot.draw(props{:})     ; title(title2) ; 
subplot(233); res_miq.draw(props{:})         ; title(title3) ; 
subplot(234); res_ioq_conn_f.draw(props{:})  ; title(title4) ; 
subplot(235); res_ioq_cot_f.draw(props{:})   ; title(title5) ; 
subplot(236); res_miq_flipped.draw(props{:}) ; title(title6) ; 

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




