%% ========================================================================
%  Plot the JL distortion:
%   1. histogram of pair-wise distortions
%   2. # of pairs with distortion > eps as epsilon increases
%  ========================================================================

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex')

OUT_FOLDER = 'results/JL_distortion';
mkdir(OUT_FOLDER);
EPS = 0.5;
JLFAC = 24;
USE_GPU = true;
SEED = 112;
SAVE = false;
LOAD = false;
FONT_SIZE = 26;

disp('Loading mesh...')
%fp = '../../../data/bunny.off';
fp = '../../../data/multires/multires/bunnyBotsch_multires/mesh_1.off';
%fp = '../../../data/bunnies/bunny_99k_faces.off';
%fp = '../../../data/ashish_nob/santa.off';
[~, meshname, ~] = fileparts(fp);
m = Mesh(fp); ne = m.nE; nv = m.nV; nf = m.nF;

if USE_GPU
    gd = gpuDevice();
else
    gd = [];
end

disp('Mesh loaded')

%  ------------------------------------------------------------------------
%  Compute exact and approx resistance distance
%  ------------------------------------------------------------------------
if ~LOAD
    disp('Approximating resistance distance...')
    tic
    Rtilde = resistance_distance(m, EPS, JLFAC, USE_GPU); 
    if ~isempty(gd), wait(gd); end; elapsed_Rtilde = toc
    Rtilde = gather(Rtilde);
    Rtilde = squareform(Rtilde);

    disp('Calculating exact resistance distance (might take a while)...')
    tic
    R = resistance_distance(m, 0, JLFAC, false, false);
    if ~isempty(gd), wait(gd); end; elapsed_R = toc
    R = gather(R);
    %fprintf('|R - Rtilde| = %.4g\n', norm(R - squareform(Rtilde), 'fro'));
else
    disp('loading workspace...')
    filename = fullfile(OUT_FOLDER, 'workspace');
    load(filename);
    disp('workspace loaded')
end

% ------------------------------------------------------------------------
% Plot histogram of pair-wise distortions
% ------------------------------------------------------------------------
%distortions = R ./ Rtilde;

% Calculate the pair-wise distortions
mask = tril(true(size(R)), -1);
Rl = R(mask);

%EPS = 1.3:-0.3:0.7;
%KS = 500:1000:2500;
%KS = nv:-1000:nv-2000;
%dims = [1000, 750, 500, 250, 100, 75, 50, 25, 10, 5, 2];
dims = [1000, 500, 250];
N = length(dims);
distortions = {};
for i = N:-1:1
    %eps = EPS(i);
    k = dims(i);
    disp(k)
    %disp(eps)
    tic
    Rtilde = resistance_distance(m, EPS, JLFAC, USE_GPU, k); 
    if ~isempty(gd), wait(gd); end; elapsed_Rtilde = toc
    Rtilde = gather(Rtilde);
    %Rtilde = squareform(Rtilde);
    
    distortions{i} = Rtilde ./ Rl;
    %mask = tril(true(size(R)), -1);
    %distortions = distortions(mask);
end

%% Plot them
bins = linspace(0.7, 1.3, 50);
%bins = 100;
colors = linspecer(N);
fh = figure; hold on
lgds = {};
hh = {};
for i = 1:N
    k = dims(i);
    disp(k)
    %subplot(1,N,i)
    [nn, xx] = hist(distortions{i}, bins);
    %nn = nn / sum(nn);
    hh{i} = bar(xx, nn, 'hist');
    set(hh{i}, 'FaceColor', colors(i, :))
    set(hh{i}, 'FaceAlpha', 0.5)
    
    %pd = fitdist(nn', 'Normal');
    %yyfit = pdf(pd, xx);
    %plot(xx, yyfit, '--k', 'LineWidth', 2);
    
    f = fit(xx(:), nn(:), 'gauss2');
    h = plot(f, '--k');
    set(h, 'LineWidth', 1.5);
    
    %hdl{i} = histfit(distortions{i});
    %hdl{i} = hist(distortions{i}, bins);
    %hdl{i}(1).FaceColor = colors(i, :);
    %hdl{i}(2).Color = 'k';
    %lgds{i} = ['\epsilon = ', num2str(eps)];
    %lgds{end+1} = ['k = ', num2str(k)];
    eps = sqrt( (24*log(nv)) / k );
    lgds{end+1} = sprintf('$k = %d, \\epsilon=%.2f$', k, eps);
    lgds{end+1} = '';
end

%%
plt = Plot(fh, true);
%plt.Title = sprintf('Histogram of pair-wise ratios R(i,j) / Rtilde(i,j)');
plt.ShowBox = 'off';
plt.FontSize = FONT_SIZE;
plt.Legend = lgds;
lh = legend;
lh.FontSize = floor(FONT_SIZE*1.5/2)-3;
%legend(lgds, 'FontSize', 12);
plt.XLim = [0.7, 1.3];
%plt.YLim = [0, 0.11]; 
plt.XLabel = '$R_{ij} \\tilde(R)_{\\epsilon}_{ij}$';
plt.YLabel = '';
plt.YTick = {};
plt.YTickLabel = {};
plt.XLabel='$R_{ij} / \tilde R_{ij}$';
hold off

set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');

if SAVE
    filename = fullfile(OUT_FOLDER, [meshname '_jl_disto_hist']);
    %print(gcf, filename, '-dpng', '-r300')
    %export_fig([filename '.pdf'])
    print(gcf, '-dpdf', [filename '.pdf']);
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

disp('Done hists')

return

%%
% ------------------------------------------------------------------------
% Plot the # of pairs with distortion > (1 +- epsilon) as epsilon increases
% ------------------------------------------------------------------------
%EPSILONS = 0.3:0.2:2;
%EPSILONS = 0.07:0.01:0.4;
%EPSILONS = 1:-0.1:0.1;
dims = [1000, 750, 500, 250, 100, 75, 50, 25, 10, 5, 2];
REPS = 30;
%DC = zeros(length(EPSILONS), REPS); % Distortion count
DC = zeros(length(dims), REPS); % Distortion count
R = gather(R);
mask = tril(true(size(R)), -1);
Rl = R(mask);
total = nchoosek(nv, 2);
tol = 0.1;
progressbar('Dimensions', 'Reps')
%for i = 1:length(EPSILONS)
    %eps = EPSILONS(i);
for i = 1:length(dims)
    k = dims(i);
    disp(eps)
    
    for j = 1:REPS
        disp('Resistance distance...')
        [Rtilde, Ztilde] = resistance_distance(m, EPS, JLFAC, USE_GPU, k);
        Rtilde = gather(Rtilde);
        
        disp('Distortion...')
        %Rtilde = squareform(Rtilde);
        distortions2 = Rtilde ./ Rl;
        %mask = tril(true(size(R)), -1);
        %distortions = distortions(mask);
        dcount = nnz(distortions2 > 1+tol | distortions2 < 1-tol);
        DC(i, j) = dcount / total;
    
        frac2 = j / REPS;
        frac1 = ((i-1) + frac2) / length(dims);
        progressbar(frac1, frac2)
    end
end

%%

fh = figure;
DCmean = mean(DC, 2);
DCstd = std(DC, 0, 2);
errorbar(dims, DCmean, DCstd, 'color', 'k')
hold on; plot(dims, DCmean, 'color', colors(1,:), 'LineWidth', 1.5); hold off
plt = Plot(fh, true);
%plt.Title = sprintf('Percentage of pairs with distortions greater than 1+-%.2g. Original dimension = 100k.', tol);
plt.ShowBox = 'off';
plt.FontSize = FONT_SIZE;
plt.XLabel = 'Reduced Dimension $k$';
plt.YLabel = '\%';
plt.YTickLabel = 100*plt.YTick;

%plt.Legend = lgds;

%xlabel('Projected dimension')
%ylabel('%')
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');

if SAVE
    filename = fullfile(OUT_FOLDER, [meshname '_jl_disto_perc']);
    %print(gcf, filename, '-dpng', '-r300')
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

%% Save workspace
if ~LOAD && SAVE
    filename = fullfile(OUT_FOLDER, [meshname '_workspace']);
    save(filename, 'R', 'DC', 'fp', 'dims', 'colors', 'FONT_SIZE', 'distortions', 'lgds', '-v7.3');
    %exclude_vars = '^(?!(m|res_tc|V|F|alpha|beta)$). ';
    %save(filename,'-regexp', exclude_vars)
end