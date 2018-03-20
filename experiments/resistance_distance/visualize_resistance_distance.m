%% ========================================================================
%  Visualize the resistance distance and the approximated resistance
%  distance on some mesh (relative to some vertex).
%  ========================================================================
clear all
plotting_defaults;
out_folder = 'results';
mkdir(out_folder);
EPS = 0.3;
JLFAC = 24;
USE_GPU = true;
SEED = 112;
SAVE = true;
FONT_SIZE = 40;
FONT = 'Computer Modern Roman';
VERT = 135;
fp = '../../../data/elephant_r.off';
%fp = '../../../data/bunnies/bunny_26k_faces.off';
%fp = '../../../data/multires/multires/bunnyBotsch_multires/mesh_2.off';
%fp = '../../../data/3holes_lo.off';
[~, meshname, ~] = fileparts(fp);
m = Mesh(fp); ne = m.nE;

if USE_GPU
    gd = gpuDevice();
else
    gd = [];
end

%%
tic
Rtilde = resistance_distance(m, EPS, JLFAC, USE_GPU); 
if ~isempty(gd), wait(gd); end; elapsed_Rtilde = toc
Rtilde = gather(Rtilde);
tic
R = resistance_distance(m, 0, JLFAC, true);
if ~isempty(gd), wait(gd); end; elapsed_R = toc
fprintf('|R - Rtilde| = %.4g\n', norm(R - squareform(Rtilde), 'fro') / norm(R, 'fro'));

% Plot
%figure
%subplot(221)
%plot_resistance(m, R, VERT); colorbar
%title({'R', sprintf('elapsed=%g', elapsed_R)})
%subplot(222)
%plot_resistance(m, squareform(Rtilde), VERT); colorbar
%title({'Rtilde', sprintf('elapsed=%g', elapsed_Rtilde)})

%%
fh = figure;
plot_resistance(m, R, VERT);
%plt.Title = {'R', sprintf('elapsed = %.4g', elapsed_R)};
plt = Plot(fh, true);
plt.Title = sprintf('R, elapsed = %.4g', elapsed_R);
plt.FontSize = FONT_SIZE;
colorbar;

%h = schemaball(R ./ abs(max(R(:))));
%schemaball(squareform(Rtilde) ./ abs(max(Rtilde(:))))

if SAVE
    filename = fullfile(out_folder, [meshname, '_R_vs_Rtilde_01']);
    export_fig([filename '.png']);
    export_fig([filename '.pdf'])
    saveas(gcf, filename, 'fig')
    %filename = fullfile(out_folder, [meshname, '_R_vs_Rtilde']);
    %print(gcf, filename, '-dpng', '-r300')
end

fh = figure;
plot_resistance(m, squareform(Rtilde), VERT);
plt = Plot(fh, true);
plt.Title = sprintf('$\\tilde R$, elapsed = %.4g, $\\epsilon=%.2g$', elapsed_Rtilde, EPS);
plt.FontSize = FONT_SIZE;
colorbar;

if SAVE
    filename = fullfile(out_folder, [meshname, '_R_vs_Rtilde_02']);
    export_fig([filename '.png']);
    export_fig([filename '.pdf'])
    saveas(gcf, filename, 'fig')
    %filename = fullfile(out_folder, [meshname, '_R_vs_Rtilde']);
    %print(gcf, filename, '-dpng', '-r300')
end

%%
close all
nv = m.nV;
Rtilde = {R};
names = {sprintf('$n=%d, T=%.2f$', nv, elapsed_R)};
i = 2;
progressbar
EPSILONS = [0.3, 0.5, 1];
max_r = max(R(:));
min_r = min(R(:));

for eps = EPSILONS
    rng(SEED)
    tic
    Rtilde{i} = resistance_distance(m, eps, JLFAC, USE_GPU); 
    if ~isempty(gd), wait(gd); end; elapsed = toc;
    Rtilde{i} = gather(Rtilde{i});
    Rt = Rtilde{i};
    max_r = max(max_r, max(Rt(:)));
    min_r = min(min_r, min(Rt(:)));
    %names{i} = {sprintf('$\\tilde R$, $\\epsilon=%g$', eps), ...
    %    sprintf('err=%.4g', norm(R-squareform(Rtilde{i}), 'fro') / norm(R, 'fro')), ...
    %    sprintf('elapsed=%.4g', elapsed)};
    k = round(24*log(nv)/eps^2);
    names{i} = sprintf('$k = %d, \\epsilon=%.1f, T = %.2f$', k, eps, elapsed);
    progressbar((i-1) / length(0.2:0.1:1))
    i = i + 1;
end
%%
close all
VIEW = [4.7324, 73.0606];
for i = 1:numel(names)
    fh = figure;
    % a hack to place a title below the mesh using legend
    plot3(m.V(1,1),m.V(1,2),m.V(1,3),'w') 
    X = Rtilde{i};
    if size(X, 2) == 1
        X = squareform(X);
    end
    plot_resistance(m, X, VERT, 'View', VIEW);
    
    plt = Plot(fh, true);
    %plt.Title = names{i};
    %plt.FontSize = FONT_SIZE;
    
    
    %lh = legend(names{i}, 'Location', 'south');
    %lh.Box = 'off';
    %lh.FontSize = FONT_SIZE;
    %%lh.Position=[0.2558 0.0409 0.2602 0.0289];
    %lh.Position = [0.2118 0.0409 0.2082 0.0289];
    
    %plt.FontName = FONT;
    %if i == 1
    %    colorbar('Location', 'West')
    %    caxis([min_r, max_r])
    %end
    colormap(linspecer)
    set(gcf, 'WindowStyle', 'docked')
    set(gcf,'color','w');
    
    %camzoom(2);
    %ax = gca;
    %ax.Position = [0.8125    0.2865    6.0000    3.0000];
    %ax.CameraPosition = [-0.3825   -2.3639    6.4491];
    %ax.CameraPosition = [0.3164, -2.5537, 6.3595];
    
    if SAVE
        name = sprintf('vis_res_dist_%s_%04d', meshname, i);
        filename = fullfile(out_folder, name);
        export_fig([filename '.png']);
        export_fig([filename '.pdf']);
        saveas(gcf, filename, 'fig')
        create_title([filename, '_title'], names{i}, FONT_SIZE);
    end
end
fh = figure;
ch = colorbar;
caxis([gather(min_r), gather(max_r)])
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
axis off
colormap(linspecer)
plt=Plot(fh, true);
plt.FontSize = 15;
if SAVE
    name = sprintf('colorbar_%s', meshname);
    filename = fullfile(out_folder, name);
    export_fig([filename '.png']);
    export_fig([filename '.pdf']);
    saveas(gcf, filename, 'fig')
end

% figure
% for i = 1:numel(names)
%     subplot(3, ceil(numel(names) / 3), i)
%     X = Rtilde{i};
%     if size(X, 2) == 1
%         X = squareform(X);
%     end
%     plot_resistance(m, X, VERT)
%     title(names{i})
% end
% 
% if SAVE
%     filename = fullfile(out_folder, [meshname, '_R_vs_Rtilde_varying_eps']);
%     print(gcf, filename, '-dpng', '-r300')
%     saveas(gcf, filename, 'fig')
%     filename = fullfile(out_folder, 'workspace');
%     save(filename);
% end