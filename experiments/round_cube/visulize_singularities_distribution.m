%% ========================================================================
%  Run JL IOQ with various epsilons and visualize the distribution of the
%  singularities (as a function on the vertices).
%  ========================================================================

%% Setup
%fp = '../../../data/bunny.off';
%fp = '../../../data/round_cuber.off';
plotting_defaults;
fp = '../../../data/rounded_cube_keenan.off';

FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0];

out_folder = fullfile('results');
mkdir(out_folder);
%if ~exist(out_folder, 'dir')
%    mkdir(out_folder);
%end
[~, filename, ~] = fileparts(fp);
filename = fullfile(out_folder, [filename, '_sing_dist']);

SAVE = false;
REPS = 300;
SAME_SEED = false;
SEED = 112;
USE_GPU = true;
FONT_SIZE=26;

m = Mesh(fp);
V = m.F; F = m.F; nv = m.nV; ne = m.nE;
d0 = get_exterior_derivatives(m);
L = d0' * d0;
Lp = inv(gpuArray(single(full(L))) + 1/nv) - 1/nv;
EPS = 0.3:0.1:1;

%%
results = {};
names = {};
n_experiments = length(EPS) + 1;
funcs = zeros(n_experiments, nv);
progressbar
prog = 0;
for r = 1:REPS
    fprintf('%d // %d\n', r, REPS);
    counter = 1;
    
    %% ------------------------------------------------------------------------
    %  True resistance
    %  ------------------------------------------------------------------------
    if SAME_SEED, rng(SEED); end
    [alpha_p1, beta_p1, elapsed_total1] = IOQ_highgenus_gpu(V, F, ...
                            'highg_method', 'genus0', ...
                            'InvMethod', 'CholMexInv', ...
                            'bsx', false, ...
                            'Mesh', m, ...
                            'UseGPU', USE_GPU, ...
                            'LaplacianPInv', Lp);
    k = [alpha_p1; beta_p1];
    funcs(counter, :) = funcs(counter, :) + k';
    counter = counter + 1;
    if r == 1
        names{end+1} = {'R', ...
            sprintf('k = %d', nv), ...
            sprintf('T = %.4g', elapsed_total1)};
    end

    %% ------------------------------------------------------------------------
    %  Approx resistance
    %  ------------------------------------------------------------------------
    for eps = EPS
        fprintf('\teps = %g\n', eps);
        if SAME_SEED, rng(SEED); end
        [alpha_p2, beta_p2, elapsed_total2] = IOQ_highgenus_gpu(V, F, ...
            'InvMethod', 'ApproxResistance', ...
            'highg_method', 'genus0', ...
            'Iterations', 30, ...
            'UseGPU', USE_GPU, ...
            'Mesh', m, ...
            'JLEps', eps, ...
            'Colamd', true);

        k = [alpha_p2; beta_p2];
        funcs(counter, :) = funcs(counter, :) + k';
        counter = counter + 1;
        if r == 1
            names{end+1} = {['Rtilde, eps=', num2str(eps)], ...
                sprintf('k = %d', round(24*log(nv)/eps^2)), ...
                sprintf('T = %.4g', elapsed_total2)};
        end
        
        prog = prog + 1;
        progressbar(prog / (REPS * length(EPS)))
    end

end

%% ------------------------------------------------------------------------
%  MIQ
%  ------------------------------------------------------------------------
%[theta, p, kmiq, R, local_frames, frame_diffs, elapsed, pFixed] = ...
%    NRoSy_mex_bin(fp, FACE0, GVEC, DEGREE);
%Emiq = E_MIQ(m, theta, frame_diffs, p, DEGREE);
tic
[m, elapsed] = nrosy_mex(fp, FACE0, GVEC, DEGREE);
elapsedMIQ = toc;
inds = m.vert_sing(:,1);
sing = m.vert_sing(:,2);
funcs(counter, inds) = sing;
counter = counter + 1;
names{end+1} = {'MIQ', sprintf('T = %.2f', elapsedMIQ)};

%% ------------------------------------------------------------------------
%  Plot results
%  ------------------------------------------------------------------------
%plotting_defaults;
%set(0,'DefaultFigureColormap',cbrewer('div','Spectral',64));
close all;
n_splots = length(names);
to_plot = [1, 2, 4, 7];
CM = cbrewer('seq', 'Blues', 256);

cam = [];
cam.pba = [1 1 1];
cam.dar = [1 1 1];
cam.cva = 10.9785;
cam.cuv = [0 0 1];
cam.ct = [-0.3612 0.0652 -0.2913];
cam.cp = [-2.6770 -2.9529 1.9051];
opt = {'PlotField', false, ...
       'PlotSing', false, ...
       'FaceAlpha', 1, ...
       'EdgeAlpha', 0, ...
       'Dock', true, ...
       'Colormap', CM, ...
       'Caxis', 'auto', ...
       'Camera', cam};
for i = to_plot
    %subplot(3, ceil(n_splots/3), i)
    fh = figure;
    plot3(m.V(1,1),m.V(1,2),m.V(1,3),'w')
    m.draw(funcs(i,:)'/sum(funcs(i,:)), opt{:})
    plt = Plot(fh, true);
    %plt.Title = names{i};
    if i == 1
        lgd = sprintf('$n = %d$', nv);
    else
        eps = EPS(i-1);
        lgd = sprintf('$k = %d, \\epsilon = %.1f$', round(24*log(nv)/eps^2), eps);
    end
    lh = legend(lgd, 'Location',  'south');
    %legend_handles{end+1} = lh;
    %lh.Position = [0.2493, 0.0359, 0.1401, 0.0289];
    lh.Box = 'off';
    lh.FontSize = FONT_SIZE*1.5;
    %lh.Position(2) = lh.Position(2)-0.015;
    plt.FontSize = FONT_SIZE*1.5;
    lh.Position = [0.2810    0.0100    0.4599    0.0828];
    %title(names{i});
    %set(gcf,'units','normalized','outerposition',[0 0 1 1]) % fullscreen so that text shows
    %colorbar
    
    set(gcf, 'WindowStyle', 'docked')
    set(gcf,'color','w');
    if SAVE
        %filename = fullfile(out_folder, 'ioq_cmp_timing_eps');
        filename = sprintf('round_cube_jl_0%d', i);
        export_fig([filename '.png']);
        export_fig([filename '.pdf']);
        %print(filename, '-dpng', '-r300')
        saveas(gcf, filename, 'fig')
        %save(filename, 'funcs', 'names', 'm')
    end

    %plt = Plot(fh, true);
    
    %plt.Title = names{i};
    %plt.FONT_SIZE = FONT_SIZE;
    %set(TitleH, 'Position', [0.5, 1], ...
  %'VerticalAlignment', 'bottom', ...
  %'HorizontalAlignment', 'center')
    %set(get(gca,'title'),'Position',[5.5 0.4 1.00011])
    %colormap(linspecer)
    
end
% if SAVE
%      %filename = fullfile(out_folder, 'ioq_cmp_timing_eps');
%      print(filename, '-dpng', '-r300')
%      saveas(gcf, filename, 'fig')
%      save(filename, 'funcs', 'names', 'm')
% end

%
figure
%m.draw
ch = colorbar;
ch.FontSize=FONT_SIZE;
F = funcs(to_plot, :);
F = F ./ sum(F, 2);
caxis([min(F(:)), max(F(:))])
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
axis off
colormap(CM)
plt=Plot(fh, true);
plt.FontSize = FONT_SIZE;
if SAVE
    filename = sprintf('colorbar');
    %filename = fullfile(out_folder, name);
    export_fig([filename '.png']);
    export_fig([filename '.pdf']);
    saveas(gcf, filename, 'fig')
end
