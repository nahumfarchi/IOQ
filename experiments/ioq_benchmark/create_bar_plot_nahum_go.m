clear all
close all

stats_file = 'results\ashish_nob_02\stats';
OUT_FOLDER = 'results\ashish_nob_02\plots';
SAVE = true;
AX_FONT_SIZE = 10;
FONT_SIZE = 13;

s = load(stats_file); s = s.stats;

%for i=2:size(s,1)-2
%    names{i-1} = s{i,1};
%    for j=2:size(s,2)
%        d = s{i,j};
%        T(i-1,j-1) = d.elapsed_total;
%        E(i-1,j-1) = d.miq_energy;
%        S(i-1,j-1) = d.n_singularities;
%    end
%end

%E = stats;
%E(2:end, 2:end) = cellfun(@get_energy, stats(2:end, 2:end), 'UniformOutput', false);
algs = s(1, 2:end);
names = s(2:end, 1);
alg_to_column = containers.Map(algs, (1:length(algs)));
T = cellfun(@get_time, s(2:end, 2:end), 'UniformOutput', true);
E = cellfun(@get_energy, s(2:end, 2:end), 'UniformOutput', true);
S = cellfun(@get_sing, s(2:end, 2:end), 'UniformOutput', true);
M = cellfun(@get_nF, s(2:end, 2:end), 'UniformOutput', true);

go_column = alg_to_column('GO');
ioq_column = alg_to_column('IOQ_conn_1a');
epsilon = 0.5;
ioq_eps_column = alg_to_column(sprintf('JL_IOQ_eps%.1f',epsilon));

valid = find(~isnan(E(:,ioq_eps_column)));
T = T(valid, :);
E = E(valid, :);
S = S(valid, :);
M = M(valid, :);
names = names(valid);

[ss,is] = sort(E(:,go_column)+T(:,go_column));
S = S(is,:);
T = T(is,:);
E = E(is,:);
names = names(is);

columns = [ioq_column, ioq_eps_column];
for i = 1:length(columns)
    j = columns(i);
    impr(:,i) = E(:,go_column).*(S(:,go_column)+1)./E(:,j)./(S(:,j)+1);
    imprT(:,i) = T(:,go_column)./T(:,j);
end

col = linspecer(3); cc = col(2,:); col(2,:) = col(3,:); col(3,:) = cc;
%%
close all;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex');
set(0,'defaulttextinterpreter', 'latex')

% E NS
[ss, is] = sort(impr(:, 1));
impr = impr(is, :);
names = names(is);
%N = round(size(S, 1) / 2);
N = round(length(names) / 2);

fh=figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3.2, 6];
plt.ShowBox = 'off';
c = categorical(names(1:N));
c = reordercats(c, names(1:N));
h = barh(c,impr(1:N, :),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
line([1 1], get(gca, 'ylim'), 'LineWidth', 1.5, 'color', 'k', 'LineStyle', '--');
%lh = legend('$(E \cdot ns)_{\textrm{MIQ}} / (E \cdot ns)_{\textrm{IOQ}}$', '$(E \cdot ns)_{\textrm{MIQ}} / (E \cdot ns)_{\textrm{IOQe}}$', 'MIQ');
%lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'GO');
lh.FontSize = FONT_SIZE;
%plt.Legend = {'$(E*ns)_M / (E*ns)_I$', '$(E*ns)_M / (E*ns)_J$'};
%plt.XLim = [0, 7];
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlim([0 3.5]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['E_NS_benchmark_01_GO']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh=figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [6, 6];
plt.ShowBox = 'off';
c = categorical(names(N+1:end));
c = reordercats(c, names(N+1:end));
h = barh(c,impr(N+1:end, :),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
line([1 1], get(gca, 'ylim'), 'LineWidth', 1.5, 'color', 'k', 'LineStyle', '--');
%lh = legend('$(E \cdot ns)_{\textrm{MIQ}} / (E \cdot ns)_{\textrm{IOQ}}$', '$(E \cdot ns)_{\textrm{MIQ}} / (E \cdot ns)_{\textrm{IOQe}}$', 'MIQ');
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'GO');
lh.FontSize = FONT_SIZE;
%plt.Legend = {'$(E*ns)_M / (E*ns)_I$', '$(E*ns)_M / (E*ns)_J$'};
%plt.XLim = [0, 7];
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlim([0 7]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['E_NS_benchmark_02_GO']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end
    
impr = impr(is, :);
names = names(is);

fprintf('Median improvement (IOQ, IOQe), (%.2f,%.2f)\n', median(impr(:,1)), median(impr(:,2)));

return

%% Timing
for column = [go_column, go_column, ioq_column, ioq_eps_column]
    fprintf('min,max timing %s (%.2f, %.2f) seconds \n', ...
        algs{column}, min(T(:,column)), max(T(:,column)));
end

[ss, is] = sort(imprT(:, 2));
imprT = imprT(is, :);
names = names(is);
N = round(length(names) / 2);

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 4];
plt.ShowBox = 'off';
c = categorical(names(1:N));
c = reordercats(c, names(1:N));
h = barh(c,imprT(1:N,:),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
%lh = legend('$T_{\textrm{MIQ}} / T_{\textrm{IOQ}}$','$T_{\textrm{MIQ}} / T_{\textrm{IOQe}}$');
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'GO');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlim([0 .1]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['T_benchmark_01_GO']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 4];
plt.ShowBox = 'off';
c = categorical(names(N+1:end));
c = reordercats(c, names(N+1:end));
h = barh(c,imprT(N+1:end,:),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
%lh = legend('$T_{\textrm{MIQ}} / T_{\textrm{IOQ}}$','$T_{\textrm{MIQ}} / T_{\textrm{IOQe}}$');
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'GO');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlim([0 .6]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['T_benchmark_02_GO']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

impr = imprT(is, :);
names = names(is);

return

%% Singularities
[ss, is] = sort(S(:, ioq_column));
S = S(is, :);
names = names(is);
%N = round(size(S, 1) / 2);
N = round(length(names) / 2);

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 4];
plt.ShowBox = 'off';
c = categorical(names(1:N));
c = reordercats(c, names(1:N));
columns = [ioq_column, ioq_eps_column, go_column];
h = barh(c,S(1:N,columns)); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
h(3).FaceColor = col(3,:);
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'GO');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Number of singularities |S|')
xlim([0 100]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['NS_benchmark_01_GO']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 4];
plt.ShowBox = 'off';
c = categorical(names(N+1:end));
c = reordercats(c, names(N+1:end));
columns = [ioq_column, ioq_eps_column, go_column];
h = barh(c,S(N+1:end,columns)); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
h(3).FaceColor = col(3,:);
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'GO');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Number of singularities |S|')
xlim([50 1000]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['NS_benchmark_02_GO']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

S = S(is, :);
names = names(is);

%% Energy
[ss, is] = sort(E(:, ioq_column));
E = E(is, :);
names = names(is);
N = round(length(names) / 2);

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 4];
plt.ShowBox = 'off';
c = categorical(names(1:N));
c = reordercats(c, names(1:N));
columns = [ioq_column, ioq_eps_column, go_column];
h = barh(c,E(1:N,columns)); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
h(3).FaceColor = col(3,:);
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'GO');
xlabel('Enregy |E|')
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlim([0 500]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['E_benchmark_01_GO']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 4];
plt.ShowBox = 'off';
c = categorical(names(N+1:end));
c = reordercats(c, names(N+1:end));
columns = [ioq_column, ioq_eps_column, go_column];
h = barh(c,E(N+1:end,columns)); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
h(3).FaceColor = col(3,:);
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'GO');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Enregy |E|')
xlim([0 500]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['E_benchmark_02_GO']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

E = E(is, :);
names = names(is);

%%

set(groot, 'defaultAxesTickLabelInterpreter', 'none'); 
set(groot, 'defaultLegendInterpreter', 'none');
set(0,'defaulttextinterpreter', 'none')



