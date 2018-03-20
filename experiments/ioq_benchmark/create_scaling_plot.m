%clear all
close all

stats_file = 'results\bunnyBotsch_multires\stats';
OUT_FOLDER = 'results\bunnyBotsch_multires\plots';
SAVE = false;
AX_FONT_SIZE = 26;
FONT_SIZE = 26;
N_ALGS = 3;

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
T  = cellfun(@get_time, s(2:end, 2:end), 'UniformOutput', true);
E  = cellfun(@get_energy, s(2:end, 2:end), 'UniformOutput', true);
S  = cellfun(@get_sing, s(2:end, 2:end), 'UniformOutput', true);
NE = cellfun(@get_ne, s(2:end, 2:end), 'UniformOutput', true);

miq_column = alg_to_column('MIQ');
ioq_column = alg_to_column('IOQ_conn_1a');
ioq_eps_column = alg_to_column('JL_IOQ_eps0.5');
go_column = alg_to_column('GO');

valid = find(~isnan(E(:,ioq_eps_column)));
T  = T(valid, :);
E  = E(valid, :);
S  = S(valid, :);
NE = NE(valid, :);
names = names(valid);

% number of faces
NF = cellfun(@(x) x.nF, s(2:end, miq_column+1));
NF = NF(valid);
[NF,is] = sort(NF);
S = S(is,:);
T = T(is,:);
E = E(is,:);
names = names(is);

% columns = [ioq_column, ioq_eps_column];
% for i = 1:length(columns)
%     j = columns(i);
%     impr(:,i) = E(:,miq_column).*(S(:,miq_column)+1)./E(:,j)./(S(:,j)+1);
%     imprT(:,i) = T(:,miq_column)./T(:,j);
% end

%col = linspecer(3); %cc = col(2,:); col(2,:) = col(3,:); col(3,:) = cc;
%%
close all;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex');
set(0,'defaulttextinterpreter', 'latex')

opt = [];
%opt.Colors = {col(1,:), col(2,:), col(3,:)};
opt.Colors = linspecer(N_ALGS);
opt.Colors = opt.Colors([2, 1, 3], :);
%opt.Colors = cbrewer('qual', 'Set1', 8);
%opt.LineWidth = [2.5, 2.5, 2.5];
opt.LineWidth = [1.5, 1.5, 1.5];
opt.LineStyle = {'-', ':', '--'};
%opt.Markers = {'o', 'x', 's'};
opt.ShowBox = 'off';
opt.FontSize = FONT_SIZE;
opt.XLabel = 'Number of faces';
opt.YLabel = 'Energy';
opt.LegendBox = 'on';
DIM = [3, 4];
opt.BoxDim = DIM;
opt.XScale = 'linear';
opt.YScale = 'linear';

fh = figure;
%E = E ./ NE;
rows = [1, 7, 8:14];
plot(NF(rows), E(rows, miq_column), ...
     NF(rows), E(rows, ioq_column), ...
     NF(rows), E(rows, ioq_eps_column));
setPlotProp(opt);
plt = Plot(fh, true);
plt.ShowBox = 'off';
plt.BoxDim = DIM;
plt.XLim = [1, 2.4096e+05];
plt.YLim = [ 0.3552   97.0413];
legend('MIQ', 'IOQ', '$IOQ_{\epsilon}$', 'Location', 'northwest')
%set(gcf, 'WindowStyle', 'docked')

if SAVE
    filename = fullfile(OUT_FOLDER, 'scaling_energy');
    %print(gcf, filename, '-dpng', '-r300')
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure;
opt.YLabel = 'Time (s)';
opt.YMinorTick = 'off';
opt.XScale = 'log';
opt.YScale = 'log';
plot(NF, T(:, miq_column), ...
     NF, T(:, ioq_column), ...
     NF, T(:, ioq_eps_column));
setPlotProp(opt);
plt = Plot(fh, true);
plt.ShowBox = 'off';
plt.BoxDim = DIM;
%xlim([5, 0.2*10^6]);
%plt.Legend = {'MIQ', 'IOQ', '$IOQ_{\epsilon}$'};
legend('MIQ', 'IOQ', '$IOQ_{\epsilon}$', 'Location', 'northwest');
set(gcf, 'WindowStyle', 'docked')
%opt.YScale = 'linear';
opt.YMinorTick = 'on';


if SAVE
    filename = fullfile(OUT_FOLDER, 'scaling_timing');
    %print(gcf, filename, '-dpng', '-r300')
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure;
opt.YLabel = 'Number of Singularities';
plot(NF, S(:, miq_column), ...
     NF, S(:, ioq_column), ...
     NF, S(:, ioq_eps_column));
setPlotProp(opt); 
plt = Plot(fh, true);
plt.ShowBox = 'off';
plt.BoxDim = DIM;
%xlim([5, 0.2*10^6]);
legend('MIQ', 'IOQ', '$IOQ_{\epsilon}$', 'Location', 'northwest')

set(gcf, 'WindowStyle', 'docked')

if SAVE
    filename = fullfile(OUT_FOLDER, 'scaling_NS');
    %print(gcf, filename, '-dpng', '-r300')
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

%%
% fh=figure; hold on;
% c = categorical(names);
% h = barh(c,impr); set(gca, 'TickLabelInterpreter', 'none');
% h(1).FaceColor = col(1,:);
% h(2).FaceColor = col(2,:);
% line([1 1], get(gca, 'ylim'), 'LineWidth', 1.5, 'color', 'k', 'LineStyle', '--');
% legend('$(E \cdot ns)_M / (E \cdot ns)_I$', '$(E \cdot ns)_M / (E \cdot ns)_J$', 'MIQ');
% %plt.Legend = {'$(E*ns)_M / (E*ns)_I$', '$(E*ns)_M / (E*ns)_J$'};
% %plt.XLim = [0, 7];
% ax = gca; set(gca,'XGrid','on','Ydir','reverse');
% ax.FontSize = AX_FONT_SIZE;
% xlim([0 7]);
% set(gcf, 'WindowStyle', 'docked')
% set(gcf,'color','w');
% if SAVE
%     filename = fullfile(OUT_FOLDER, ['E_NS_benchmark']);
%     export_fig([filename '.pdf'])
%     export_fig([filename '.png'])
%     saveas(gcf, filename, 'fig')
% end
% 
% figure; hold on;
% c = categorical(names);
% h = barh(c,imprT); set(gca, 'TickLabelInterpreter', 'none');
% h(1).FaceColor = col(1,:);
% h(2).FaceColor = col(2,:);
% legend('$T_M / T_I$','$T_M / T_J$');
% ax = gca; set(gca,'XGrid','on','Ydir','reverse');
% ax.FontSize = AX_FONT_SIZE;
% xlim([0 1.2]);
% set(gcf, 'WindowStyle', 'docked')
% set(gcf,'color','w');
% if SAVE
%     filename = fullfile(OUT_FOLDER, ['T_benchmark']);
%     export_fig([filename '.pdf'])
%     export_fig([filename '.png'])
%     saveas(gcf, filename, 'fig')
% end
% 
% figure; hold on;
% c = categorical(names);
% columns = [ioq_column, ioq_eps_column, miq_column];
% h = barh(c,S(:,columns)); set(gca, 'TickLabelInterpreter', 'none');
% h(1).FaceColor = col(1,:);
% h(2).FaceColor = col(2,:);
% h(3).FaceColor = col(3,:);
% legend('\#s IOQ', '\#s $IOQ_{\epsilon}$','\#s MIQ');
% ax = gca; set(gca,'XGrid','on','Ydir','reverse');
% ax.FontSize = AX_FONT_SIZE;
% xlim([0 50]);
% set(gcf, 'WindowStyle', 'docked')
% set(gcf,'color','w');
% if SAVE
%     filename = fullfile(OUT_FOLDER, ['NS_benchmark']);
%     export_fig([filename '.pdf'])
%     export_fig([filename '.png'])
%     saveas(gcf, filename, 'fig')
% end
% 
% figure; hold on;
% c = categorical(names);
% columns = [ioq_column, ioq_eps_column, miq_column];
% h = barh(c,E(:,columns)); set(gca, 'TickLabelInterpreter', 'none');
% h(1).FaceColor = col(1,:);
% h(2).FaceColor = col(2,:);
% h(3).FaceColor = col(3,:);
% legend('E IOQ', 'E $IOQ_{\epsilon}$','E MIQ');
% ax = gca; set(gca,'XGrid','on','Ydir','reverse');
% ax.FontSize = AX_FONT_SIZE;
% xlim([0 500]);
% set(gcf, 'WindowStyle', 'docked')
% set(gcf,'color','w');
% if SAVE
%     filename = fullfile(OUT_FOLDER, ['E_benchmark']);
%     export_fig([filename '.pdf'])
%     export_fig([filename '.png'])
%     saveas(gcf, filename, 'fig')
% end

%%
close all;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex');
set(0,'defaulttextinterpreter', 'latex')

colors = linspecer(N_ALGS);
colors = colors([2, 1, 3], :);

% opt.LineWidth = [1.5, 1.5, 1.5];
% opt.LineStyle = {'-', ':', '--'};
% opt.Markers = {'o', 'x', 's'};
% opt.ShowBox = 'off';
% opt.FontSize = FONT_SIZE;
% opt.XLabel = 'Number of faces';
% opt.YLabel = 'Energy';
% opt.LegendBox = 'on';
% opt.BoxDim = [6, 3];
% opt.XScale = 'linear';
% opt.YScale = 'linear';

fh = figure;
%E = E ./ NE;
%plt.Colors = colors;
%plt.ShowBox = 'off';
plot(NF, E(:, miq_column), ...
     NF, E(:, ioq_column), ...
     NF, E(:, ioq_eps_column));
legend('MIQ', 'IOQ', '$IOQ_{\epsilon}$')

plt = Plot();
set(gcf, 'WindowStyle', 'docked')



%%

if SAVE
    filename = fullfile(OUT_FOLDER, 'scaling_energy');
    %print(gcf, filename, '-dpng', '-r300')
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure;
plt = Plot(fh, true);
plt.YLabel = 'Elapsed Time';
plt.YMinorTick = 'off';
plt.XScale = 'log';
plt.YScale = 'log';
plot(NF, T(:, miq_column), ...
     NF, T(:, ioq_column), ...
     NF, T(:, ioq_eps_column));
legend('MIQ', 'IOQ', '$IOQ_{\epsilon}$', 'Location', 'northwest')
%setPlotProp(opt);
set(gcf, 'WindowStyle', 'docked')
%opt.YScale = 'linear';
%opt.YMinorTick = 'on';

if SAVE
    filename = fullfile(OUT_FOLDER, 'scaling_timing');
    %print(gcf, filename, '-dpng', '-r300')
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure;
plt = Plot(fh, true);
plt.YLabel = 'Number of Singularities';
plot(NF, S(:, miq_column), ...
     NF, S(:, ioq_column), ...
     NF, S(:, ioq_eps_column));
legend('MIQ', 'IOQ', '$IOQ_{\epsilon}$')
%setPlotProp(opt);
set(gcf, 'WindowStyle', 'docked')

if SAVE
    filename = fullfile(OUT_FOLDER, 'scaling_NS');
    %print(gcf, filename, '-dpng', '-r300')
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

%%
close all;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex'); 
set(groot, 'defaultLegendInterpreter', 'latex');
set(0,'defaulttextinterpreter', 'latex')
colors = linspecer(3);
DIM = [6, 3];

fh = figure;
%rows = [1, 7, 8:14];
ph = plot(NF(rows), E(rows, miq_column), ...
     NF(rows), E(rows, ioq_column), ...
     NF(rows), E(rows, ioq_eps_column));
ph(1).Color = colors(2,:);
ph(2).Color = colors(1,:);
ph(3).Color = colors(3,:);
plt = Plot(fh, true);
plt.XLabel = 'Number of Faces';
plt.YLabel = 'Energy';
plt.ShowBox = 'off';
plt.BoxDim = DIM;
plt.Colors = colors;
plt.Legend = {'MIQ', 'IOQ', 'IOQe'};
plt.LegendBox = 'on';
plt.LegendLoc = 'northwest';
plt.YLim = [0 45];
plt.FontSize = 18;
plt.XScale = 'log';
set(gcf, 'WindowStyle', 'docked')
if SAVE
    filename = fullfile(OUT_FOLDER, 'scaling_energy');
    %print(gcf, filename, '-dpng', '-r300')
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure;
ph = plot(NF, T(:, miq_column), ...
     NF, T(:, ioq_column), ...
     NF, T(:, ioq_eps_column), '--');
ph(1).Color = colors(2,:);
ph(2).Color = colors(1,:);
ph(3).Color = colors(3,:);
plt = Plot(fh, true);
plt.XLabel = 'Number of Faces';
plt.YLabel = 'Timing (s)';
plt.ShowBox = 'off';
plt.BoxDim = DIM;
plt.Colors = colors;
plt.Legend = {'MIQ', 'IOQ', 'IOQe'};
plt.LegendBox = 'on';
plt.LegendLoc = 'northwest';
plt.YLim = [0 45];
plt.FontSize = 18;
plt.XScale = 'log';
plt.YScale = 'log';
set(gcf, 'WindowStyle', 'docked')

%%

set(groot, 'defaultAxesTickLabelInterpreter', 'none'); 
set(groot, 'defaultLegendInterpreter', 'none');
set(0,'defaulttextinterpreter', 'none')

%%
close all
COLORS = linspecer(3);
COLORS = {COLORS(2,:), COLORS(1, :), COLORS(3, :)};
BOX_DIM = [3, 4];

fh = figure;
cols = [1:5, 7:length(NF)];
plot(NF(cols), E(cols, miq_column), ...
     NF(cols), E(cols, ioq_column), ...
     NF(cols), E(cols, ioq_eps_column));
plt = Plot(fh, false);
plt.Colors = COLORS;
plt.Legend = {'MIQ', 'IOQ', 'IOQe'};
plt.YScale = 'log';
plt.BoxDim = BOX_DIM;
plt.LineStyle = {'-', '-', ':'};
plt.ShowBox = 'off';
plt.LegendBox = 'on';
plt.XLabel = 'Number of Faces';
plt.YLabel = 'Energy';
plt.XLim = [0, 1.1129e5];
plt.LegendLoc = 'northwest';
plt.YLim = [0,1.6e2];
set(gcf, 'WindowStyle', 'docked')
filename = fullfile(OUT_FOLDER, 'scaling_energy.pdf');
%plt.export(filename);
if SAVE, export_fig(filename); end

fh = figure;
plot(NF, T(:, miq_column), ...
     NF, T(:, ioq_column), ...
     NF, T(:, ioq_eps_column));
plt = Plot(fh, false);
plt.Colors = COLORS;
plt.Legend = {'MIQ', 'IOQ', 'IOQe'};
plt.XScale = 'log';
plt.YScale = 'log';
plt.BoxDim = BOX_DIM;
plt.LineStyle = {'-', '-', ':'};
plt.ShowBox = 'off';
plt.LegendBox = 'on';
plt.XLabel = 'Number of Faces';
plt.YLabel = 'Time (s)';
plt.XLim = [0, 2e5];
plt.XTick = 10.^(2:5)
plt.LegendLoc = 'northwest';
set(gcf, 'WindowStyle', 'docked')
filename = fullfile(OUT_FOLDER, 'scaling_timing.pdf');
%plt.export(filename);
if SAVE, export_fig(filename); end


fh = figure;
plot(NF, S(:, miq_column), ...
     NF, S(:, ioq_column), ...
     NF, S(:, ioq_eps_column));
plt = Plot(fh, false);
plt.Colors = COLORS;
plt.Legend = {'MIQ', 'IOQ', 'IOQe'};
plt.XScale = 'log';
%plt.YScale = 'log';
plt.BoxDim = BOX_DIM;
plt.LineStyle = {'-', '-', ':'};
plt.ShowBox = 'off';
plt.LegendBox = 'on';
plt.XLabel = 'Number of Faces';
plt.YLabel = 'Number of Singularities';
plt.XLim = [0, 2e5];
plt.XTick = 10.^(2:5)
plt.LegendLoc = 'northwest';
set(gcf, 'WindowStyle', 'docked')
filename = fullfile(OUT_FOLDER, 'scaling_NS.pdf');
%plt.export(filename);
if SAVE, export_fig(filename); end

 


