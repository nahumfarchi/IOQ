clear all
close all

stats_file = 'results\ashish_nob_02\stats';
OUT_FOLDER = 'results\ashish_nob_02\plots';
SAVE = false;
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

miq_column = alg_to_column('MIQ');
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
onames = names;
oT = T;
oE = E;
oS = S;
oM = M;

[ss,is] = sort(E(:,miq_column)+T(:,miq_column));
S = S(is,:);
T = T(is,:);
E = E(is,:);
names = names(is);

columns = [ioq_column, ioq_eps_column];
for i = 1:length(columns)
    j = columns(i);
    impr(:,i) = E(:,miq_column).*(S(:,miq_column)+1)./E(:,j)./(S(:,j)+1);
    imprT(:,i) = T(:,miq_column)./T(:,j);
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
%lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'MIQ');
lh.FontSize = FONT_SIZE;
%plt.Legend = {'$(E*ns)_M / (E*ns)_I$', '$(E*ns)_M / (E*ns)_J$'};
%plt.XLim = [0, 7];
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlim([0 3.5]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['E_NS_benchmark_01']);
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
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'MIQ');
lh.FontSize = FONT_SIZE;
%plt.Legend = {'$(E*ns)_M / (E*ns)_I$', '$(E*ns)_M / (E*ns)_J$'};
%plt.XLim = [0, 7];
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlim([0 7]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['E_NS_benchmark_02']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end
    
impr = impr(is, :);
names = names(is);

fprintf('Median improvement (IOQ, IOQe), (%.2f,%.2f)\n', median(impr(:,1)), median(impr(:,2)));

close all

names = onames;
T = oT;
E = oE;
S = oS;
M = oM;


%% Timing
for column = [miq_column, go_column, ioq_column, ioq_eps_column]
    fprintf('min,max,median timing %s (%.2f, %.2f, %.2f) seconds \n', ...
        algs{column}, min(T(:,column)), max(T(:,column)), median(T(:,column)));
end

% % [ss, is] = sort(imprT(:, 2));
% % imprT = imprT(is, :);
% % names = names(is);
% % N = round(length(names) / 2);
% % 
% % fh = figure; hold on;
% % plt = Plot(fh, true);
% % plt.BoxDim = [3, 4];
% % plt.ShowBox = 'off';
% % c = categorical(names(1:N));
% % c = reordercats(c, names(1:N));
% % h = barh(c,imprT(1:N,:),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
% % h(1).FaceColor = col(1,:);
% % h(2).FaceColor = col(2,:);
% % %lh = legend('$T_{\textrm{MIQ}} / T_{\textrm{IOQ}}$','$T_{\textrm{MIQ}} / T_{\textrm{IOQe}}$');
% % lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'MIQ');
% % lh.FontSize = FONT_SIZE;
% % ax = gca; set(gca,'XGrid','on','Ydir','reverse');
% % ax.FontSize = AX_FONT_SIZE;
% % xlim([0 .1]);
% % set(gcf, 'WindowStyle', 'docked')
% % set(gcf,'color','w');
% % if SAVE
% %     filename = fullfile(OUT_FOLDER, ['T_benchmark_01']);
% %     export_fig([filename '.pdf'])
% %     export_fig([filename '.png'])
% %     saveas(gcf, filename, 'fig')
% % end
% % 
% % fh = figure; hold on;
% % plt = Plot(fh, true);
% % plt.BoxDim = [3, 4];
% % plt.ShowBox = 'off';
% % c = categorical(names(N+1:end));
% % c = reordercats(c, names(N+1:end));
% % h = barh(c,imprT(N+1:end,:),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
% % h(1).FaceColor = col(1,:);
% % h(2).FaceColor = col(2,:);
% % %lh = legend('$T_{\textrm{MIQ}} / T_{\textrm{IOQ}}$','$T_{\textrm{MIQ}} / T_{\textrm{IOQe}}$');
% % lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'MIQ');
% % lh.FontSize = FONT_SIZE;
% % ax = gca; set(gca,'XGrid','on','Ydir','reverse');
% % ax.FontSize = AX_FONT_SIZE;
% % xlim([0 .6]);
% % set(gcf, 'WindowStyle', 'docked')
% % set(gcf,'color','w');
% % if SAVE
% %     filename = fullfile(OUT_FOLDER, ['T_benchmark_02']);
% %     export_fig([filename '.pdf'])
% %     export_fig([filename '.png'])
% %     saveas(gcf, filename, 'fig')
% % end
% 
% impr = imprT(is, :);
% names = names(is);
% 
% return


%% Timing IOQ
names = onames;
T = oT;
E = oE;
S = oS;
M = oM;
[ss, is] = sort(T(:, ioq_column));
T = T(is, :);
names = names(is);
%N = round(size(S, 1) / 2);
N = 37; %round(length(names) / 2);

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 4.5];
plt.ShowBox = 'off';
c = categorical(names(1:N));
c = reordercats(c, names(1:N));
columns = [ioq_column];
h = barh(c,T(1:N,columns),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
lh = legend('IOQ');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Timing in seconds')
xlim([0 30]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['T_benchmark_01_IOQ']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 4.5];
plt.ShowBox = 'off';
N2 = 54;
c = categorical(names(N+1:N2));
c = reordercats(c, names(N+1:N2));
columns = [ioq_column];
h = barh(c,T(N+1:N2,columns),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
lh = legend('IOQ');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Timing in seconds')
xlim([160 610]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['T_benchmark_02_IOQ']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 4.5];
plt.ShowBox = 'off';
c = categorical(names(N2+1:end));
c = reordercats(c, names(N2+1:end));
columns = [ioq_column];
h = barh(c,T(N2+1:end,columns),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
lh = legend('IOQ');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Timing in seconds')
xlim([900 1800]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['T_benchmark_03_IOQ']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end


%% Timing IOQe
names = onames;
T = oT;
E = oE;
S = oS;
M = oM;
[ss, is] = sort(T(:, ioq_eps_column));
T = T(is, :);
names = names(is);
%N = round(size(S, 1) / 2);
N = 42; %round(length(names) / 2);

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 5];
plt.ShowBox = 'off';
c = categorical(names(1:N));
c = reordercats(c, names(1:N));
columns = [ioq_eps_column];
h = barh(c,T(1:N,columns),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(2,:);
lh = legend(['IOQe, ' num2str(epsilon)]);
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Timing in seconds')
xlim([0 16]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['T_benchmark_01_IOQe']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 5];
plt.ShowBox = 'off';
N2 = 82;
c = categorical(names(N+1:N2));
c = reordercats(c, names(N+1:N2));
columns = [ioq_eps_column];
h = barh(c,T(N+1:N2,columns),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(2,:);
lh = legend(['IOQe, ' num2str(epsilon)]);
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Timing in seconds')
xlim([15 105]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['T_benchmark_02_IOQe']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 5];
plt.ShowBox = 'off';
c = categorical(names(N2+1:end));
c = reordercats(c, names(N2+1:end));
columns = [ioq_eps_column];
h = barh(c,T(N2+1:end,columns),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(2,:);
lh = legend(['IOQe, ' num2str(epsilon)]);
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Timing in seconds')
xlim([100 710]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['T_benchmark_03_IOQe']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end



T = T(is, :);
names = names(is);


%% Singularities
names = onames;
T = oT;
E = oE;
S = oS;
M = oM;

[ss, is] = sort(S(:, ioq_column));
S = S(is, :);
names = names(is);
%N = round(size(S, 1) / 2);
N = 40;%round(length(names) / 2);

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 6];
plt.ShowBox = 'off';
c = categorical(names(1:N));
c = reordercats(c, names(1:N));
columns = [ioq_column, ioq_eps_column, miq_column];
h = barh(c,S(1:N,columns),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
h(3).FaceColor = col(3,:);
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'MIQ');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Number of singularities')
xlim([0 95]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['NS_benchmark_01']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

N2 = 75;
fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [6, 6];
plt.ShowBox = 'off';
c = categorical(names(N+1:N2));
c = reordercats(c, names(N+1:N2));
columns = [ioq_column, ioq_eps_column, miq_column];
h = barh(c,S(N+1:N2,columns),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
h(3).FaceColor = col(3,:);
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'MIQ');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Number of singularities')
xlim([50 280]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['NS_benchmark_02']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 6];
plt.ShowBox = 'off';
c = categorical(names(N2+1:end));
c = reordercats(c, names(N2+1:end));
columns = [ioq_column, ioq_eps_column, miq_column];
h = barh(c,S(N2+1:end,columns),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
h(3).FaceColor = col(3,:);
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'MIQ');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Number of singularities')
xlim([100 1000]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['NS_benchmark_03']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

S = S(is, :);
names = names(is);

names = onames;
T = oT;
E = oE;
S = oS;
M = oM;

close all

SAVE = true;

%% Energy
[ss, is] = sort(E(:, ioq_column));
E = E(is, :);
names = names(is);
N = round(length(names) / 3);

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 6];
plt.ShowBox = 'off';
c = categorical(names(1:N));
c = reordercats(c, names(1:N));
columns = [ioq_column, ioq_eps_column, miq_column];
h = barh(c,E(1:N,columns),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
h(3).FaceColor = col(3,:);
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'MIQ');
xlabel('Enregy')
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlim([0 45]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['E_benchmark_01']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [3, 6];
plt.ShowBox = 'off';
N2 = 2*N;
c = categorical(names(N+1:N2));
c = reordercats(c, names(N+1:N2));
columns = [ioq_column, ioq_eps_column, miq_column];
h = barh(c,E(N+1:N2,columns),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
h(3).FaceColor = col(3,:);
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'MIQ');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Enregy')
xlim([20 90]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['E_benchmark_02']);
    export_fig([filename '.pdf'])
    export_fig([filename '.png'])
    saveas(gcf, filename, 'fig')
end

fh = figure; hold on;
plt = Plot(fh, true);
plt.BoxDim = [6, 6];
plt.ShowBox = 'off';
c = categorical(names(N2+1:end));
c = reordercats(c, names(N2+1:end));
columns = [ioq_column, ioq_eps_column, miq_column];
h = barh(c,E(N2+1:end,columns),'edgecolor','none'); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
h(3).FaceColor = col(3,:);
lh = legend('IOQ', ['IOQe, ' num2str(epsilon)], 'MIQ');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
xlabel('Enregy')
xlim([50 500]);
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
if SAVE
    filename = fullfile(OUT_FOLDER, ['E_benchmark_03']);
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



