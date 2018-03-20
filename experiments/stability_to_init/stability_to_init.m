%% ========================================================================
%  Stability to initialization.
%  ========================================================================
plotting_defaults;
FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = true;
SAVE = false;
LOAD = true;
FONT_SIZE = 26;
%OUT_FOLDER = 'results/laurent_hand_nob_r';
OUT_FOLDER = 'results/bunnyBotsch/';
mkdir(OUT_FOLDER);

%fp = '../../../data/rounded_cube_keenan_r_3k.off';
%fp = '../../../data/bunny.off';
%fp = '../../../data/multires/multires/dragon_multires/mesh_10.off';
fp = '../../../data/multires/multires/bunnyBotsch_multires/mesh_4.off';
%fp = '../../../data/multires/multires/gearbox_multires/mesh_4.off';
%fp = '../../../data/laurent_hand_nob_r.off';
[~, meshname, ~] = fileparts(fp);
%OUT_FOLDER = fullfile(OUT_FOLDER, meshname);
%mkdir(OUT_FOLDER);
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;

d0 = get_exterior_derivatives(m);
L = d0'*d0;
try
    Lp = inv(gpuArray(single(full(L))) + 1/nv) - 1/nv;
catch
    Lp = invChol_mex(single(full(L)) + 1/nv) - 1/nv;
end

% run ioq
[alpha, beta] = IOQ_highgenus_gpu(...
    V, F, ...
    'highg_method', 'option1a', ...
    'Mesh', m, ...
    'UseGPU', USE_GPU, ...
    'LaplacianPInv', Lp);
k = [alpha; beta];
res_ioq = TCODS(m, ...
    'k', k, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);

% run miq
[res_miq, elapsed_miq] = ...
    nrosy_mex(fp, FACE0, GVEC, DEGREE);



REPS = 30;
N_INIT_SING = 0:10:100;

Sxx = zeros(length(N_INIT_SING), 1);
S = zeros(length(N_INIT_SING), REPS);
E = zeros(length(N_INIT_SING), REPS);
VFUNCSi = zeros(length(N_INIT_SING), nv); % initial singularities distribution
VFUNCSf = zeros(length(N_INIT_SING), nv); % final singularities distribution

if ~LOAD
    progressbar('Singularities', 'Reps')
    for i = 1:length(N_INIT_SING)
        ns = N_INIT_SING(i);
        disp(i)
        for j = 1:REPS
            [alpha, beta, ~, ~, out] = IOQ_highgenus_gpu(...
                V, F, ...
                'highg_method', 'option1a', ...
                'Mesh', m, ...
                'UseGPU', USE_GPU, ...
                'LaplacianPInv', Lp, ...
                'NSingularities', ns);

            k = [alpha; beta];
            res_tc = TCODS(m, ...
                'k', k, ...
                'f0', FACE0, ...
                'theta0', THETA0, ...
                'degree', DEGREE, ...
                'CreateFField', true, ...
                'Duplicate', true, ...
                'gConstraintVec', GVEC);

            S(i, j) = nnz(alpha);
            E(i, j) = res_tc.miq_energy;
            VFUNCSi(i, :) = VFUNCSi(i, :) + out.init_alpha';
            VFUNCSf(i, :) = VFUNCSf(i, :) + alpha';

            frac2 = j / REPS;
            frac1 = ((i-1) + frac2) / length(N_INIT_SING);
            progressbar(frac1, frac2)
        end
        Sxx(i) = out.n_init_sing;
    end
    progressbar(1, 1)
else
    filename = fullfile(OUT_FOLDER, 'workspace');
    load(filename)
end

if ~LOAD && SAVE
    filename = fullfile(OUT_FOLDER, 'workspace');
    save(filename, 'S', 'Sxx', 'E', 'VFUNCSi', 'VFUNCSf', 'fp');
    %exclude_vars = '^(?!(m|res_tc|V|F|alpha|beta)$). ';
    %save(filename,'-regexp', exclude_vars)
end

%% Plot number of final singularities against number of initial singularities
opt = [];
%opt.Colors = {col(1,:), col(2,:), col(3,:)};
opt.Colors = linspecer(2);
%opt.Colors = cbrewer('qual', 'Set1', 8);
%opt.LineWidth = [2.5, 2.5, 2.5];
opt.LineWidth = [2, 2, 2];
opt.LineStyle = {'-', '--', '--'};
opt.Markers = {'', '', 's'};
opt.ShowBox = 'off';
opt.FontSize = FONT_SIZE;
opt.XLabel = 'Initial Singularities';
opt.YLabel = 'Final Singularities';
opt.LegendBox = 'on';
%opt.BoxDim = [12, 6];
opt.YLim = [22, 45];

fh = figure; hold on
Smean = mean(S, 2);
Sstd = std(S, 0, 2);
errorbar(Sxx, Smean, Sstd)
plot(Sxx, res_miq.n_vert_sing * ones(length(Sxx), 1));
%color = linspecer(1);
%h = barwitherr(Sstd, Smean, 'FaceColor', color);
%bar(Sstd, Smean);
%title('Number of final singularities against number of initial singularities')
%xlabel('# initial singularities')
%ylabel('# final singularities')
legend('IOQ', 'MIQ')
set(gca,'XTickLabel', Sxx)

setPlotProp(opt);
set(gcf, 'WindowStyle', 'docked')

%opt = [];
%opt.XLabel = '# of initial singularities';
%opt.YLabel = '# of final singularities';
%opt.Title = 'Number of final singularities against number of initial singularities';
%%opt.FileName = fullfile(OUT_FOLDER, 'final_sing_against_init_sing.png';
%setPlotProp(opt, fh);

% plt = Plot(fh, true);
% plt.XLabel = 'Initial # of cones';
% plt.YLabel = 'Final # of cones';
% %plt.Title = 'Number of final singularities against number of initial singularities';
% plt.XMinorTick = 'off';
% plt.YMinorTick = 'off';
% plt.ShowBox = 'off';
% plt.FontSize = FONT_SIZE;
% plt.Colors = linspecer(2);
% set(gcf, 'WindowStyle', 'docked')

if SAVE
    filename = fullfile(OUT_FOLDER, 'final_sing_against_init_sing');
    %print(gcf, filename, '-dpng', '-r300')
    %plt.export([filename '.pdf']);
    export_fig([filename '.pdf']);
    export_fig([filename '.png']);
    saveas(gcf, filename, 'fig')
end

%% Plot final energy against number of initial singularities
fh = figure; hold on
Emean = mean(E, 2);
Estd = std(E, 0, 2);
errorbar(Sxx, Emean, Estd);
plot(Sxx, res_miq.miq_energy * ones(length(Sxx), 1));

opt.YLim = ylim;
opt.YLabel = 'Final Energy';
lh = legend('IOQ', 'MIQ');
%lh.Position(2)=0.31;
set(gca,'XTickLabel', Sxx)

setPlotProp(opt);
set(gcf, 'WindowStyle', 'docked')
lh.Position = [0.5566    0.3100    0.1283    0.1094];

% color = linspecer(1);
% %errorbar(Sxx, Emean, Estd)
% h = barwitherr(Estd, Emean, 'FaceColor', color);
% %xlabel('# initial singularities')
% %ylabel('Energy')
% plt = Plot(fh, true);
% plt.XLabel = 'Initial # of cones';
% plt.YLabel = 'Final energy';
% %plt.Title = 'Number of final singularities against number of initial singularities';
% plt.XMinorTick = 'off';
% plt.YMinorTick = 'off';
% plt.ShowBox = 'off';
% plt.FontSize = FONT_SIZE;
% plt.XTickLabel = Sxx;

if SAVE
    filename = fullfile(OUT_FOLDER, 'energy_against_init_sing');
    %print(gcf, filename, '-dpng', '-r300')
    %saveas(gcf, filename, 'fig')
    export_fig([filename '.pdf']);
    export_fig([filename '.png']);
    saveas(gcf, filename, 'fig')
end

%% Heatmap plot of the distribution of final singularities
m = Mesh(fp);
% %figure
% 
% figure
% set(0,'DefaultFigureColormap',cbrewer('div','Spectral',64));
% [ha, pos] = tight_subplot(1, 2, 0, 0, 0);
% 
% colormap(linspecer)
% row = 10;
% 
% %subplot(121)
% axes(ha(1));
% fi = VFUNCSi(row, :)';
% fi = fi / sum(fi);
% m.draw('Func', fi, 'FaceAlpha', 1, 'EdgeAlpha', 0, 'PlotField', false, 'PlotSing', false)
% colorbar
% title(['Initial singularities distribution, init\_ns=', num2str(Sxx(row))])
% 
% %subplot(122)
% axes(ha(2));
% ff = VFUNCSf(row, :)';
% ff = ff / sum(ff);
% m.draw('Func', ff, 'FaceAlpha', 1, 'EdgeAlpha', 0, 'PlotField', false, 'PlotSing', false)
% colorbar
% title('Final singualrities distribution')
% 
% if SAVE
%     filename = fullfile(OUT_FOLDER, 'sing_dist');
%     %print(gcf, filename, '-dpng', '-r300')
%     %saveas(gcf, filename, 'fig')
%     export_fig([filename '.pdf']);
%     export_fig([filename '.png']);
%     saveas(gcf, filename, 'fig')
% end

%% Heatmap plot of the distribution of final singularities with wfigs
%m.draw(m.V(:,1))
% cam = [];
% cam.pba = [1.3876 1.0944 1];
% cam.dar = [1 1 1];
% cam.cva = 6.0386;
% cam.cuv = [0.0292 0.9996 0.0030];
% cam.ct = [-0.0174 0.1009 0.0012];
% cam.cp = [1.4289 0.0583 0.1363];

% cam = [];
% cam.pba = [1.3876 1.0944 1];
% cam.dar = [1 1 1];
% cam.cva = 6.0386;
% cam.cuv = [-0.0179 0.9998 0.0039];
% cam.ct = [-0.0174 0.1009 0.0012];
% cam.cp = [-1.4692 0.0747 0.0601];

cam = [];
cam.pba = [1.387598 1.094399 1.000000 ];
cam.dar = [1 1 1 ];
cam.cva = [6.038600 ];
cam.cuv = [-0.017885 0.999832 0.003900 ];
cam.ct = [-0.017400 0.100900 0.001200 ];
cam.cp = [-0.905408 0.080529 1.151376 ];

close all;

%
m = Mesh(fp);
row = 10;
%VIEW = [-359.2861, 89.4902];
%VIEW = [159.9307, 34.2773];
VIEW = [87.9139, 5.4902];
filename = fullfile(OUT_FOLDER, 'sing_dist2');
%CM = linspecer(256);
%CM = CM(floor(size(CM,1)/2):end, :);
%CM = CM(1:floor(size(CM,1)/2), :);
CM = cbrewer('div', 'RdBu', 256);
%CM = cbrewer('div', 'Spectral', 256);
opt = {'PlotField', false, ...
       'PlotSing', false, ...
       'FaceAlpha', 1, ...
       'EdgeAlpha', 0, ...
       'Dock', true, ...
       'View', VIEW, ...
       'Colormap', CM, ...
       'Caxis', 'auto', ...
       'Camera', cam};
%f1 = VFUNCSi(row,:)'; f1 = f1 - mean(f1); f1 = f1 / sum(abs(f1));
%f2 = VFUNCSf(row,:)'; f2 = f2 - mean(f2); f2 = f2 / sum(abs(f2));
f1 = VFUNCSi(row,:)'; f2 = VFUNCSf(row,:)';
%m.montage(filename, [f1, f2], 1, opt{:});
figure; m.draw('Func', f1, opt{:}); th = title('(a)', 'FontSize', FONT_SIZE*1.5); th.Position=[0,0.017,0];
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
export_fig(sprintf('%s_%04d.pdf', filename, 1));
figure; m.draw('Func', f2, opt{:}); th = title('(b)', 'FontSize', FONT_SIZE*1.5); th.Position=[0,0.017,0];
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
export_fig(sprintf('%s_%04d.pdf', filename, 2));

%
%
% ioq
%close all
fh = figure; hold on
% a hack to place a title below the mesh using legend
plot3(m.V(1,1),m.V(1,2),m.V(1,3),'w') 
%opt = {'PlotField', false, 'PlotSing', true, 'Dock', false, 'View', [0.8, 90]};
%opt = {'PlotField', false, 'PlotSing', true, 'Dock', false, 'View', [-359.2861, 89.4902]};
opt = {'PlotField', false, ...
       'PlotSing', true, ...
       'FaceAlpha', 1, ...
       'EdgeAlpha', 0, ...
       'Dock', true, ...
       'View', VIEW, ...
       'Colormap', CM, ...
       'Caxis', 'auto', ...
       'Camera', cam};
res_ioq.draw(opt{:});
title_ioq = {'(c) IOQ', sprintf('$E = %.2f, |S| = %d$', res_ioq.miq_energy, res_ioq.n_vert_sing)};
th = title(title_ioq, 'FontSize', FONT_SIZE*1.5);
th.Position=[0,0.002,0];
%th.Position = [0.3000    0.008   -0.0012];

%lgd = legend(title_ioq, 'Location', 'south');
%lgd.Box = 'off';
%lgd.FontSize = 35;
%lgd.Position=[0.3176 0.25 0.3999 0.1253];

%text(5, 0.4, title_ioq);
%title(title_ioq)
%plt = Plot(fh, true);
%plt.Title = title_ioq;
%plt.FontSize = FONT_SIZE;
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, 'ioq');
    export_fig([filename '.pdf']);
    export_fig([filename '.png']);
    saveas(gcf, filename, 'fig');
end

% miq
fh = figure; hold on
plot3(m.V(1,1),m.V(1,2),m.V(1,3),'w')
res_miq.draw(opt{:});
%title_miq = sprintf('MIQ, $E = %.2f, ns = %d$', res_miq.miq_energy, res_miq.n_vert_sing);
title_miq = {'(d) MIQ', sprintf('$E = %.2f, |S| = %d$', res_miq.miq_energy, res_miq.n_vert_sing)};
%text(0, 0, title_miq);
th = title(title_miq, 'FontSize', FONT_SIZE*1.5);
th.Position=[0,0.002,0];
%th.Position = [0.3000    0.008   -0.0012];

%lgd = legend(title_miq, 'Location', 'south');
%lgd.Box = 'off';
%lgd.FontSize = 35;
%lgd.Position=[0.3176 0.25 0.3999 0.1253];

%plt = Plot(fh, true);
%plt.Title = title_miq;
%plt.FontSize = FONT_SIZE;
%plt.LegendBox = 'off';
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
%
if SAVE
    filename = fullfile(OUT_FOLDER, 'miq');
    export_fig([filename '.pdf']);
    export_fig([filename '.png']);
    saveas(gcf, filename, 'fig');
end

%% close up
cam = [];
cam.pba = [1.387598 1.094399 1.000000 ];
cam.dar = [1 1 1 ];
cam.cva = [6.038600 ];
cam.cuv = [-0.017437 0.999832 -0.005682 ];
cam.ct = [-0.016780 0.054684 0.005432 ];
cam.cp = [0.701610 0.067403 0.039028 ];

close all;

%
m = Mesh(fp);
row = 10;
%VIEW = [-359.2861, 89.4902];
%VIEW = [159.9307, 34.2773];
VIEW = [87.9139, 5.4902];
filename = fullfile(OUT_FOLDER, 'sing_dist2_closeup');
%CM = linspecer(256);
%CM = CM(floor(size(CM,1)/2):end, :);
%CM = CM(1:floor(size(CM,1)/2), :);
CM = cbrewer('div', 'RdBu', 256);
%CM = cbrewer('div', 'Spectral', 256);
opt = {'PlotField', false, ...
       'PlotSing', false, ...
       'FaceAlpha', 1, ...
       'EdgeAlpha', 0, ...
       'Dock', true, ...
       'View', VIEW, ...
       'Colormap', CM, ...
       'Caxis', 'auto', ...
       'Camera', cam, ...
       'MarkerSize', 160};
%f1 = VFUNCSi(row,:)'; f1 = f1 - mean(f1); f1 = f1 / sum(abs(f1));
%f2 = VFUNCSf(row,:)'; f2 = f2 - mean(f2); f2 = f2 / sum(abs(f2));
f1 = VFUNCSi(row,:)'; f2 = VFUNCSf(row,:)';
%m.montage(filename, [f1, f2], 1, opt{:});
figure; m.draw('Func', f1, opt{:}); %th = title('(a)', 'FontSize', FONT_SIZE*1.5); th.Position=[0,0.017,0];
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
export_fig(sprintf('%s_%04d.pdf', filename, 1));
figure; m.draw('Func', f2, opt{:}); %th = title('(b)', 'FontSize', FONT_SIZE*1.5); th.Position=[0,0.017,0];
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
export_fig(sprintf('%s_%04d.pdf', filename, 2));

%%
%
% ioq
%close all
fh = figure; hold on
% a hack to place a title below the mesh using legend
plot3(m.V(1,1),m.V(1,2),m.V(1,3),'w') 
%opt = {'PlotField', false, 'PlotSing', true, 'Dock', false, 'View', [0.8, 90]};
%opt = {'PlotField', false, 'PlotSing', true, 'Dock', false, 'View', [-359.2861, 89.4902]};
opt = {'PlotField', false, ...
       'PlotSing', true, ...
       'FaceAlpha', 1, ...
       'EdgeAlpha', 0, ...
       'Dock', true, ...
       'View', VIEW, ...
       'Colormap', CM, ...
       'Caxis', 'auto', ...
       'Camera', cam, ...
       'MarkerSize', 160};
res_ioq.draw(opt{:});
%title_ioq = {'(c) IOQ', sprintf('$E = %.2f, |S| = %d$', res_ioq.miq_energy, res_ioq.n_vert_sing)};
%th = title(title_ioq, 'FontSize', FONT_SIZE*1.5);
%th.Position=[0,0.002,0];
%th.Position = [0.3000    0.008   -0.0012];

%lgd = legend(title_ioq, 'Location', 'south');
%lgd.Box = 'off';
%lgd.FontSize = 35;
%lgd.Position=[0.3176 0.25 0.3999 0.1253];

%text(5, 0.4, title_ioq);
%title(title_ioq)
%plt = Plot(fh, true);
%plt.Title = title_ioq;
%plt.FontSize = FONT_SIZE;
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, 'ioq_closeup');
    export_fig([filename '.pdf']);
    export_fig([filename '.png']);
    saveas(gcf, filename, 'fig');
end

% miq
fh = figure; hold on
plot3(m.V(1,1),m.V(1,2),m.V(1,3),'w')
res_miq.draw(opt{:});
%title_miq = sprintf('MIQ, $E = %.2f, ns = %d$', res_miq.miq_energy, res_miq.n_vert_sing);
%title_miq = {'(d) MIQ', sprintf('$E = %.2f, |S| = %d$', res_miq.miq_energy, res_miq.n_vert_sing)};
%text(0, 0, title_miq);
%th = title(title_miq, 'FontSize', FONT_SIZE*1.5);
%th.Position=[0,0.002,0];
%th.Position = [0.3000    0.008   -0.0012];

%lgd = legend(title_miq, 'Location', 'south');
%lgd.Box = 'off';
%lgd.FontSize = 35;
%lgd.Position=[0.3176 0.25 0.3999 0.1253];

%plt = Plot(fh, true);
%plt.Title = title_miq;
%plt.FontSize = FONT_SIZE;
%plt.LegendBox = 'off';
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
%
if SAVE
    filename = fullfile(OUT_FOLDER, 'miq_closeup');
    export_fig([filename '.pdf']);
    export_fig([filename '.png']);
    saveas(gcf, filename, 'fig');
end