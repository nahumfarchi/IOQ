%% ========================================================================
%  Energy per edge plots
%  x axis - edge energy
%  y axis - % of edges with energy < x
%  ========================================================================
ioq_benchmark_setup;
plotting_defaults(30, 30);
print_header(experiment_name);

filepaths = get_filepaths(out_folder_ffields, '.ffield');
n_files = numel(filepaths);

inds = cellfun(@(x) strfind(x, 'MIQ'), filepaths, 'UniformOutput', false);
inds = cellfun(@(x) ~isempty(x), inds);
filepaths_MIQ = filepaths(inds);

inds = cellfun(@(x) strfind(x, 'conn_genus0'), filepaths, 'UniformOutput', false);
inds = cellfun(@(x) ~isempty(x), inds);
filepaths_IOQ_conn_genus0 = filepaths(inds);

inds = cellfun(@(x) strfind(x, 'cot_genus0'), filepaths, 'UniformOutput', false);
inds = cellfun(@(x) ~isempty(x), inds);
filepaths_IOQ_cot_genus0 = filepaths(inds);

inds = cellfun(@(x) strfind(x, 'conn_opt1a'), filepaths, 'UniformOutput', false);
inds = cellfun(@(x) ~isempty(x), inds);
filepaths_IOQ_conn_opt1a = filepaths(inds);

inds = cellfun(@(x) strfind(x, 'conn_opt2z'), filepaths, 'UniformOutput', false);
inds = cellfun(@(x) ~isempty(x), inds);
filepaths_IOQ_conn_opt2z = filepaths(inds);

inds = cellfun(@(x) strfind(x, 'JL_IOQ_conn_eps0.5'), filepaths, 'UniformOutput', false);
inds = cellfun(@(x) ~isempty(x), inds);
filepaths_JL_IOQ_conn_eps05 = filepaths(inds);

SAVE_PLOTS = false;


%  ------------------------------------------------------------------------
%  Genus 0
%  ------------------------------------------------------------------------

% Calculate per edge energy
% MIQ
E_edge_miq_all_g0 = [];
E_miq_all_g0 = [];
progressbar
for i = 1:numel(filepaths_MIQ)
    fprintf('%d // %d\n', i, numel(filepaths_MIQ));
    progressbar(i / numel(filepaths_MIQ))
    fp = filepaths_MIQ{i};
    m = Mesh(fp);
    if m.genus > 0
        continue
    end
    E = per_edge_energy(m);
    E_edge_miq_all_g0 = [E_edge_miq_all_g0; E];
    E_miq_all_g0 = [E_miq_all_g0; m.miq_energy];
end
assert(abs(sum(E_edge_miq_all_g0) - sum(E_miq_all_g0)) < 1e-7)
fprintf('\n');

% IOQ conn
E_edge_ioq_conn_all_g0 = [];
E_ioq_conn_all_g0 = [];
progressbar
for i = 1:numel(filepaths_IOQ_conn_genus0)
    fprintf('%d // %d\n', i, numel(filepaths_IOQ_conn_genus0));
    progressbar(i / numel(filepaths_IOQ_conn_genus0))
    fp = filepaths_IOQ_conn_genus0{i};
    m = Mesh(fp);
    if m.genus > 0
        continue
    end
    E = per_edge_energy(m);
    E_edge_ioq_conn_all_g0 = [E_edge_ioq_conn_all_g0; E];
    E_ioq_conn_all_g0 = [E_ioq_conn_all_g0; m.miq_energy];
end
assert(abs(sum(E_edge_ioq_conn_all_g0) - sum(E_ioq_conn_all_g0)) < 1e-7)
fprintf('\n');

% IOQ cot
E_edge_ioq_cot_all_g0 = [];
E_ioq_cot_all_g0 = [];
progressbar
for i = 1:numel(filepaths_IOQ_cot_genus0)
    fprintf('%d // %d\n', i, numel(filepaths_IOQ_cot_genus0));
    progressbar(i / numel(filepaths_IOQ_cot_genus0))
    fp = filepaths_IOQ_cot_genus0{i};
    m = Mesh(fp);
    if m.genus > 0
        continue
    end
    E = per_edge_energy(m);
    E_edge_ioq_cot_all_g0 = [E_edge_ioq_cot_all_g0; E];
    E_ioq_cot_all_g0 = [E_ioq_cot_all_g0; m.miq_energy];
end
assert(abs(sum(E_edge_ioq_cot_all_g0) - sum(E_ioq_cot_all_g0)) < 1e-7)
fprintf('\n');

% Plot
% percentage of edges with energy < x
figure
bins = 1000;
lgd = {};

[hh1, xx1] = hist(E_edge_miq_all_g0, bins);
yy1 = cumsum(hh1) / length(E_edge_miq_all_g0) * 100;
lgd{end+1} = 'MIQ';

[hh2, xx2] = hist(E_edge_ioq_conn_all_g0, bins);
yy2 = cumsum(hh2) / length(E_edge_ioq_conn_all_g0) * 100;
lgd{end+1} = 'IOQ conn';

[hh3, xx3] = hist(E_edge_ioq_cot_all_g0, bins);
yy3 = cumsum(hh3) / length(E_edge_ioq_cot_all_g0) * 100;
lgd{end+1} = 'IOQ cot';

plot(xx1, yy1, xx2, yy2, xx3, yy3)
title('Genus 0')
xlabel('Edge energy')
ylabel('% of edges < x')
legend(lgd)
if SAVE_PLOTS
    filename = fullfile(out_folder_plots, 'per_edge_energy_genus0');
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

%% percentage of meshes with energy < x
figure
bins = 100;
lgd = {};

[hh1, xx1] = hist(E_miq_all_g0, bins);
yy1 = cumsum(hh1) / length(E_miq_all_g0) * 100;
lgd{end+1} = 'MIQ';

[hh2, xx2] = hist(E_ioq_conn_all_g0, bins);
yy2 = cumsum(hh2) / length(E_ioq_conn_all_g0) * 100;
lgd{end+1} = 'IOQ conn';

[hh3, xx3] = hist(E_ioq_cot_all_g0, bins);
yy3 = cumsum(hh3) / length(E_ioq_cot_all_g0) * 100;
lgd{end+1} = 'IOQ cot';

plot(xx1, yy1, xx2, yy2, xx3, yy3)
title('Genus 0')
xlabel('Mesh energy')
ylabel('% of meshes < x')
legend(lgd)
if SAVE_PLOTS
    filename = fullfile(out_folder_plots, 'mesh_perc_energy_genus0');
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end
%%


%  ------------------------------------------------------------------------
%  Genus > 0
%  ------------------------------------------------------------------------

% Calculate per edge energy
% MIQ
E_edge_miq_all_hg = [];
E_miq_all_hg = [];
progressbar
for i = 1:numel(filepaths_MIQ)
    fprintf('%d // %d\n', i, numel(filepaths_MIQ));
    progressbar(i / numel(filepaths_MIQ))
    fp = filepaths_MIQ{i};
    m = Mesh(fp);
    if m.genus == 0
        continue
    end
    E = per_edge_energy(m);
    E_edge_miq_all_hg = [E_edge_miq_all_hg; E];
    E_miq_all_hg = [E_miq_all_hg; m.miq_energy];
end
assert(abs(sum(E_edge_miq_all_hg) - sum(E_miq_all_hg)) < 1e-7)
fprintf('\n');

% IOQ conn 1a
E_edge_ioq_conn_all_hg_opt1a = [];
E_ioq_conn_all_hg_opt1a = [];
progressbar
for i = 1:numel(filepaths_IOQ_conn_opt1a)
    fprintf('%d // %d\n', i, numel(filepaths_IOQ_conn_opt1a));
    progressbar(i / numel(filepaths_IOQ_conn_opt1a))
    fp = filepaths_IOQ_conn_opt1a{i};
    m = Mesh(fp);
    if m.genus == 0
        continue
    end
    E = per_edge_energy(m);
    E_edge_ioq_conn_all_hg_opt1a = [E_edge_ioq_conn_all_hg_opt1a; E];
    E_ioq_conn_all_hg_opt1a = [E_ioq_conn_all_hg_opt1a; m.miq_energy];
end
assert(abs(sum(E_edge_ioq_conn_all_hg_opt1a) - sum(E_ioq_conn_all_hg_opt1a)) < 1e-7)
fprintf('\n');

% IOQ conn 2z
E_edge_ioq_conn_all_hg_opt2z = [];
E_ioq_conn_all_hg_opt2z = [];
progressbar
for i = 1:numel(filepaths_IOQ_conn_opt2z)
    fprintf('%d // %d\n', i, numel(filepaths_IOQ_conn_opt2z));
    progressbar(i / numel(filepaths_IOQ_conn_opt2z))
    fp = filepaths_IOQ_cot_genus0{i};
    m = Mesh(fp);
    if m.genus == 0
        continue
    end
    E = per_edge_energy(m);
    E_edge_ioq_conn_all_hg_opt2z = [E_edge_ioq_conn_all_hg_opt2z; E];
    E_ioq_conn_all_hg_opt2z = [E_ioq_conn_all_hg_opt2z; m.miq_energy];
end
assert(abs(sum(E_edge_ioq_conn_all_hg_opt2z) - sum(E_ioq_conn_all_hg_opt2z)) < 1e-7)
fprintf('\n');

%% Plot
figure
bins = 1000;
lgd = {};

[hh1, xx1] = hist(E_edge_miq_all_hg, bins);
yy1 = cumsum(hh1) / length(E_edge_miq_all_hg) * 100;
lgd{end+1} = 'MIQ';

[hh2, xx2] = hist(E_edge_ioq_conn_all_hg_opt1a, bins);
yy2 = cumsum(hh2) / length(E_edge_ioq_conn_all_hg_opt1a) * 100;
lgd{end+1} = 'IOQ conn 1a';

[hh3, xx3] = hist(E_edge_ioq_conn_all_hg_opt2z, bins);
yy3 = cumsum(hh3) / length(E_edge_ioq_conn_all_hg_opt2z) * 100;
lgd{end+1} = 'IOQ conn 2z';

plot(xx1, yy1, xx2, yy2, xx3, yy3)
title('Genus > 0')
xlabel('Edge energy')
ylabel('% of edges < x')
legend(lgd)
if SAVE_PLOTS
    filename = fullfile(out_folder_plots, 'per_edge_energy_highg');
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

figure
bins = 100;
lgd = {};

[hh1, xx1] = hist(E_miq_all_hg, bins);
yy1 = cumsum(hh1) / length(E_miq_all_hg) * 100;
lgd{end+1} = 'MIQ';

[hh2, xx2] = hist(E_ioq_conn_all_hg_opt1a, bins);
yy2 = cumsum(hh2) / length(E_ioq_conn_all_hg_opt1a) * 100;
lgd{end+1} = 'IOQ conn 1a';

[hh3, xx3] = hist(E_ioq_conn_all_hg_opt2z, bins);
yy3 = cumsum(hh3) / length(E_ioq_conn_all_hg_opt2z) * 100;
lgd{end+1} = 'IOQ conn 2z';

plot(xx1, yy1, xx2, yy2, xx3, yy3)
title('Genus > 0')
xlabel('Mesh energy')
ylabel('% of meshes < x')
legend(lgd)
if SAVE_PLOTS
    filename = fullfile(out_folder_plots, 'mesh_perc_energy_highg');
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

%% ========================================================================
%  energy / time plots
%  x axis - mesh name
%  y axs - energy / time
%  ========================================================================
ioq_benchmark_setup;
plotting_defaults(30, 30);
print_header(experiment_name);

stats_fp = fullfile(base_folder, 'stats.mat');
if exist(stats_fp, 'file')
    stats = load(stats_fp);
    stats = stats.stats;
else
    error('stats file does not exist')
end

exp_names = stats(1, 2:end);
name_to_col = containers.Map(exp_names, (1:length(exp_names))+1);

% sort by number of faces
miq_stats = stats(2:end, MIQ_COL);
xx = cellfun(@(x) x.nF, miq_stats);
[xx, inds] = sort(xx);
stats_sorted = stats(2:end, 1:end);
stats_sorted = stats_sorted(inds, :);
stats(2:end, 1:end) = stats_sorted;

%get_energy = @(x) x.miq_energy;
%get_time = @(x) x.elapsed_total;

is_valid = @(x) ~isempty(x) && ...
                 (~isfield(x, 'failed') || ...
                 ~x.failed);
             
is_genus0 = @(x) x.genus == 0;

SAVE_PLOTS = false;
X_FONT_SIZE = 14;

% -------------------------------------------------------------------------
% genus 0
% plot energy
% -------------------------------------------------------------------------
figure; hold on

lgds = {};
yy = stats(2:end, name_to_col('MIQ'));
inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
xx1 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
yy1 = cellfun(@(x) get_energy(x), yy);
yy1(~inds) = nan;
plot(1:length(xx1), yy1, 'X')
xticks(1:length(xx1))
xticklabels(xx1)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('MIQ'));

yy = stats(2:end, name_to_col('IOQ_conn_genus0'));
inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
xx2 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
yy2 = cellfun(@(x) get_energy(x), yy);
yy2(~inds) = nan;
plot(1:length(xx2), yy2, 'X')
xticks(1:length(xx2))
xticklabels(xx2)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('IOQ_conn_genus0'));

yy = stats(2:end, name_to_col('IOQ_cot_genus0'));
inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
xx3 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
yy3 = cellfun(@(x) get_energy(x), yy);
yy3(~inds) = nan;
plot(1:length(xx3), yy3, 'X')
xticks(1:length(xx3))
xticklabels(xx3)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('IOQ_cot_genus0'));

yy = stats(2:end, name_to_col('JL_IOQ_conn_eps0.5'));
inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
xx3 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
yy3 = cellfun(@(x) get_energy(x), yy);
yy3(~inds) = nan;
plot(1:length(xx3), yy3, 'X')
xticks(1:length(xx3))
xticklabels(xx3)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('JL_IOQ_conn_eps0.5'));

legend(lgds, 'Location', 'northwest', 'Interpreter', 'none')
title('Genus 0 Energy')

hold off

if SAVE_PLOTS
    filename = fullfile(out_folder_plots, 'energy_genus0');
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

% -------------------------------------------------------------------------
% genus 0
% plot timing
% -------------------------------------------------------------------------
figure; hold on

lgds = {};
yy = stats(2:end, name_to_col('MIQ'));
inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
xx1 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
yy1 = cellfun(@(x) get_time(x), yy);
yy1(~inds) = nan;
plot(1:length(xx1), yy1, 'X')
xticks(1:length(xx1))
xticklabels(xx1)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('MIQ'));

yy = stats(2:end, name_to_col('IOQ_conn_genus0'));
inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
xx2 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
yy2 = cellfun(@(x) get_time(x), yy);
yy2(~inds) = nan;
plot(1:length(xx2), yy2, 'X')
xticks(1:length(xx2))
xticklabels(xx2)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('IOQ_conn_genus0'));

yy = stats(2:end, name_to_col('IOQ_cot_genus0'));
inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
xx3 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
yy3 = cellfun(@(x) get_time(x), yy);
yy3(~inds) = nan;
plot(1:length(xx3), yy3, 'X')
xticks(1:length(xx3))
xticklabels(xx3)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('IOQ_cot_genus0'));

yy = stats(2:end, name_to_col('JL_IOQ_conn_eps0.5'));
inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
xx3 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
yy3 = cellfun(@(x) get_time(x), yy);
yy3(~inds) = nan;
plot(1:length(xx3), yy3, 'X')
xticks(1:length(xx3))
xticklabels(xx3)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('JL_IOQ_conn_eps0.5'));

legend(lgds, 'Location', 'northwest', 'Interpreter', 'none')
title('Genus 0 Timing')

hold off

if SAVE_PLOTS
    filename = fullfile(out_folder_plots, 'time_genus0');
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

% -------------------------------------------------------------------------
% genus > 0
% Plot energy
% -------------------------------------------------------------------------
figure; hold on

lgds = {};
yy = stats(2:end, name_to_col('MIQ'));
inds = cellfun(@(x) is_valid(x) && ~is_genus0(x), yy);
xx1 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
yy1 = cellfun(@(x) get_energy(x), yy);
yy1(~inds) = nan;
plot(1:length(xx1), yy1, 'X')
xticks(1:length(xx1))
xticklabels(xx1)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('MIQ'));

yy = stats(2:end, name_to_col('IOQ_conn_opt1a'));
inds = cellfun(@(x) is_valid(x) && ~is_genus0(x), yy);
xx2 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
%yy2 = cellfun(get_energy, yy);
yy2 = nan(size(xx2));
yy2(inds) = cellfun(@(x) get_energy(x), yy(inds));
%yy2(~inds) = nan;
plot(1:length(xx2), yy2, 'X')
xticks(1:length(xx2))
xticklabels(xx2)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('IOQ_conn_opt1a'));

yy = stats(2:end, name_to_col('IOQ_conn_opt2z'));
inds = cellfun(@(x) is_valid(x) && ~is_genus0(x), yy);
xx3 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
yy3 = nan(size(xx3));
yy3(inds) = cellfun(@(x) get_energy(x), yy(inds));
%yy3 = cellfun(get_energy, yy);
%yy3(~inds) = nan;
plot(1:length(xx3), yy3, 'X')
xticks(1:length(xx3))
xticklabels(xx3)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('IOQ_conn_opt2z'));

legend(lgds, 'Location', 'northwest', 'Interpreter', 'none')
title('High Genus Energy')

hold off

if SAVE_PLOTS
    filename = fullfile(out_folder_plots, 'energy_highg');
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

% -------------------------------------------------------------------------
% genus > 0
% plot timing
% -------------------------------------------------------------------------
figure; hold on

lgds = {};
yy = stats(2:end, name_to_col('MIQ'));
inds = cellfun(@(x) is_valid(x) && ~is_genus0(x), yy);
xx1 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
yy1 = nan(size(xx1));
yy1(inds) = cellfun(@(x) get_time(x), yy(inds));
plot(1:length(xx1), yy1, 'X')
xticks(1:length(xx1))
xticklabels(xx1)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('MIQ'));

yy = stats(2:end, name_to_col('IOQ_conn_opt1a'));
inds = cellfun(@(x) is_valid(x) && ~is_genus0(x), yy);
xx2 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
yy2 = nan(size(xx2));
yy2(inds) = cellfun(@(x) get_time(x), yy(inds));
plot(1:length(xx2), yy2, 'X')
xticks(1:length(xx2))
xticklabels(xx2)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('IOQ_conn_opt1a'));

yy = stats(2:end, name_to_col('IOQ_conn_opt2z'));
inds = cellfun(@(x) is_valid(x) && ~is_genus0(x), yy);
xx3 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
%xx1 = xx1(inds);
yy3 = nan(size(xx3));
yy3(inds) = cellfun(@(x) get_time(x), yy(inds));
plot(1:length(xx3), yy3, 'X')
xticks(1:length(xx3))
xticklabels(xx3)
xtickangle(90)
xt = get(gca, 'XTick');
set(gca, 'FontSize', X_FONT_SIZE)
lgds(end+1) = stats(1, name_to_col('IOQ_conn_opt2z'));

legend(lgds, 'Location', 'northwest', 'Interpreter', 'none')
title('High Genus Timing')

hold off

if SAVE_PLOTS
    filename = fullfile(out_folder_plots, 'time_highg');
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end











% %% ========================================================================
% %  energy / time plots
% %  x axis - mesh name
% %  y axs - energy / time
% %  ========================================================================
% ioq_benchmark_setup;
% plotting_defaults(30, 30);
% print_header(experiment_name);
% 
% stats_fp = fullfile(base_folder, 'stats.mat');
% if exist(stats_fp, 'file')
%     stats = load(stats_fp);
%     stats = stats.stats;
% else
%     error('stats file does not exist')
% end
% 
% exp_names = stats(1, 2:end);
% name_to_col = containers.Map(exp_names, (1:length(exp_names))+1);
% 
% % sort by number of faces
% miq_stats = stats(2:end, MIQ_COL);
% xx = cellfun(@(x) x.nF, miq_stats);
% [xx, inds] = sort(xx);
% stats_sorted = stats(2:end, 1:end);
% stats_sorted = stats_sorted(inds, :);
% stats(2:end, 1:end) = stats_sorted;
% 
% %get_energy = @(x) x.miq_energy;
% %get_time = @(x) x.elapsed_total;
% 
% is_valid = @(x) ~isempty(x) && ...
%                  (~isfield(x, 'failed') || ...
%                  ~x.failed);
%              
% is_genus0 = @(x) x.genus == 0;
% 
% SAVE_PLOTS = true;
% X_FONT_SIZE = 14;
% 
% % -------------------------------------------------------------------------
% % genus 0
% % plot energy
% % -------------------------------------------------------------------------
% figure; hold on
% 
% lgds = {};
% yy = stats(2:end, name_to_col('MIQ'));
% inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
% xx1 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
% %xx1 = xx1(inds);
% yy1 = cellfun(@(x) get_energy(x), yy);
% yy1(~inds) = nan;
% plot(1:length(xx1), yy1, 'X')
% xticks(1:length(xx1))
% xticklabels(xx1)
% xtickangle(90)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', X_FONT_SIZE)
% lgds(end+1) = stats(1, name_to_col('MIQ'));
% 
% yy = stats(2:end, name_to_col('IOQ_conn_genus0'));
% inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
% xx2 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
% %xx1 = xx1(inds);
% yy2 = cellfun(@(x) get_energy(x), yy);
% yy2(~inds) = nan;
% plot(1:length(xx2), yy2, 'X')
% xticks(1:length(xx2))
% xticklabels(xx2)
% xtickangle(90)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', X_FONT_SIZE)
% lgds(end+1) = stats(1, name_to_col('IOQ_conn_genus0'));
% 
% yy = stats(2:end, name_to_col('IOQ_cot_genus0'));
% inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
% xx3 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
% %xx1 = xx1(inds);
% yy3 = cellfun(@(x) get_energy(x), yy);
% yy3(~inds) = nan;
% plot(1:length(xx3), yy3, 'X')
% xticks(1:length(xx3))
% xticklabels(xx3)
% xtickangle(90)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', X_FONT_SIZE)
% lgds(end+1) = stats(1, name_to_col('IOQ_cot_genus0'));
% 
% legend(lgds, 'Location', 'northwest', 'Interpreter', 'none')
% title('Genus 0 Energy')
% 
% hold off
% 
% if SAVE_PLOTS
%     filename = fullfile(out_folder_plots, 'energy_genus0');
%     print(gcf, filename, '-dpng', '-r300')
%     saveas(gcf, filename, 'fig')
% end
% 
% % -------------------------------------------------------------------------
% % genus 0
% % plot timing
% % -------------------------------------------------------------------------
% figure; hold on
% 
% lgds = {};
% yy = stats(2:end, name_to_col('MIQ'));
% inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
% xx1 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
% %xx1 = xx1(inds);
% yy1 = cellfun(@(x) get_time(x), yy);
% yy1(~inds) = nan;
% plot(1:length(xx1), yy1, 'X')
% xticks(1:length(xx1))
% xticklabels(xx1)
% xtickangle(90)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', X_FONT_SIZE)
% lgds(end+1) = stats(1, name_to_col('MIQ'));
% 
% yy = stats(2:end, name_to_col('IOQ_conn_genus0'));
% inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
% xx2 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
% %xx1 = xx1(inds);
% yy2 = cellfun(@(x) get_time(x), yy);
% yy2(~inds) = nan;
% plot(1:length(xx2), yy2, 'X')
% xticks(1:length(xx2))
% xticklabels(xx2)
% xtickangle(90)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', X_FONT_SIZE)
% lgds(end+1) = stats(1, name_to_col('IOQ_conn_genus0'));
% 
% yy = stats(2:end, name_to_col('IOQ_cot_genus0'));
% inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
% xx3 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
% %xx1 = xx1(inds);
% yy3 = cellfun(@(x) get_time(x), yy);
% yy3(~inds) = nan;
% plot(1:length(xx3), yy3, 'X')
% xticks(1:length(xx3))
% xticklabels(xx3)
% xtickangle(90)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', X_FONT_SIZE)
% lgds(end+1) = stats(1, name_to_col('IOQ_cot_genus0'));
% 
% legend(lgds, 'Location', 'northwest', 'Interpreter', 'none')
% title('Genus 0 Timing')
% 
% hold off
% 
% if SAVE_PLOTS
%     filename = fullfile(out_folder_plots, 'time_genus0');
%     print(gcf, filename, '-dpng', '-r300')
%     saveas(gcf, filename, 'fig')
% end
% 
% % -------------------------------------------------------------------------
% % genus > 0
% % Plot energy
% % -------------------------------------------------------------------------
% figure; hold on
% 
% lgds = {};
% yy = stats(2:end, name_to_col('MIQ'));
% inds = cellfun(@(x) is_valid(x) && ~is_genus0(x), yy);
% xx1 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
% %xx1 = xx1(inds);
% yy1 = cellfun(@(x) get_energy(x), yy);
% yy1(~inds) = nan;
% plot(1:length(xx1), yy1, 'X')
% xticks(1:length(xx1))
% xticklabels(xx1)
% xtickangle(90)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', X_FONT_SIZE)
% lgds(end+1) = stats(1, name_to_col('MIQ'));
% 
% yy = stats(2:end, name_to_col('IOQ_conn_opt1a'));
% inds = cellfun(@(x) is_valid(x) && ~is_genus0(x), yy);
% xx2 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
% %xx1 = xx1(inds);
% %yy2 = cellfun(get_energy, yy);
% yy2 = nan(size(xx2));
% yy2(inds) = cellfun(@(x) get_energy(x), yy(inds));
% %yy2(~inds) = nan;
% plot(1:length(xx2), yy2, 'X')
% xticks(1:length(xx2))
% xticklabels(xx2)
% xtickangle(90)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', X_FONT_SIZE)
% lgds(end+1) = stats(1, name_to_col('IOQ_conn_opt1a'));
% 
% yy = stats(2:end, name_to_col('IOQ_conn_opt2z'));
% inds = cellfun(@(x) is_valid(x) && ~is_genus0(x), yy);
% xx3 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
% %xx1 = xx1(inds);
% yy3 = nan(size(xx3));
% yy3(inds) = cellfun(@(x) get_energy(x), yy(inds));
% %yy3 = cellfun(get_energy, yy);
% %yy3(~inds) = nan;
% plot(1:length(xx3), yy3, 'X')
% xticks(1:length(xx3))
% xticklabels(xx3)
% xtickangle(90)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', X_FONT_SIZE)
% lgds(end+1) = stats(1, name_to_col('IOQ_conn_opt2z'));
% 
% legend(lgds, 'Location', 'northwest', 'Interpreter', 'none')
% title('High Genus Energy')
% 
% hold off
% 
% if SAVE_PLOTS
%     filename = fullfile(out_folder_plots, 'energy_highg');
%     print(gcf, filename, '-dpng', '-r300')
%     saveas(gcf, filename, 'fig')
% end
% 
% % -------------------------------------------------------------------------
% % genus > 0
% % plot timing
% % -------------------------------------------------------------------------
% figure; hold on
% 
% lgds = {};
% yy = stats(2:end, name_to_col('MIQ'));
% inds = cellfun(@(x) is_valid(x) && ~is_genus0(x), yy);
% xx1 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
% %xx1 = xx1(inds);
% yy1 = nan(size(xx1));
% yy1(inds) = cellfun(@(x) get_time(x), yy(inds));
% plot(1:length(xx1), yy1, 'X')
% xticks(1:length(xx1))
% xticklabels(xx1)
% xtickangle(90)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', X_FONT_SIZE)
% lgds(end+1) = stats(1, name_to_col('MIQ'));
% 
% yy = stats(2:end, name_to_col('IOQ_conn_opt1a'));
% inds = cellfun(@(x) is_valid(x) && ~is_genus0(x), yy);
% xx2 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
% %xx1 = xx1(inds);
% yy2 = nan(size(xx2));
% yy2(inds) = cellfun(@(x) get_time(x), yy(inds));
% plot(1:length(xx2), yy2, 'X')
% xticks(1:length(xx2))
% xticklabels(xx2)
% xtickangle(90)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', X_FONT_SIZE)
% lgds(end+1) = stats(1, name_to_col('IOQ_conn_opt1a'));
% 
% yy = stats(2:end, name_to_col('IOQ_conn_opt2z'));
% inds = cellfun(@(x) is_valid(x) && ~is_genus0(x), yy);
% xx3 = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
% %xx1 = xx1(inds);
% yy3 = nan(size(xx3));
% yy3(inds) = cellfun(@(x) get_time(x), yy(inds));
% plot(1:length(xx3), yy3, 'X')
% xticks(1:length(xx3))
% xticklabels(xx3)
% xtickangle(90)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', X_FONT_SIZE)
% lgds(end+1) = stats(1, name_to_col('IOQ_conn_opt2z'));
% 
% legend(lgds, 'Location', 'northwest', 'Interpreter', 'none')
% title('High Genus Timing')
% 
% hold off
% 
% if SAVE_PLOTS
%     filename = fullfile(out_folder_plots, 'time_highg');
%     print(gcf, filename, '-dpng', '-r300')
%     saveas(gcf, filename, 'fig')
% end










%% ========================================================================
%  energy / time plots 
%  x axis - number of faces
%  y axis - energy / time
%  ========================================================================
% ioq_benchmark_setup;
% plotting_defaults(30, 30);
% stats_fp = fullfile(base_folder, 'stats.mat');
% if exist(stats_fp, 'file')
%     stats = load(stats_fp);
%     stats = stats.stats;
% else
%     error('stats file does not exist')
% end
% 
% % x axis (number of faces)
% miq_stats = stats(2:end, MIQ_COL);
% xx = cellfun(@(x) x.nF, miq_stats);
% [xx, inds] = sort(xx);
% stats_sorted = stats(2:end, 1:end);
% stats_sorted = stats_sorted(inds, :);
% stats(2:end, 1:end) = stats_sorted;
% 
% % y axis function
% get_y_field = @(x) x.elapsed_total;
% 
% is_valid = @(x) ~isempty(x) && ...
%                  (~isfield(x, 'failed') || ...
%                  ~x.failed);
%              
% is_genus0 = @(x) x.genus == 0;
% 
% MIQ_COL = 2;
% IOQ_CONN_G0_COL = 3;
% IOQ_COT_G0_COL = 4;
% OPT1A_CONN_COL = 5;
% 
% % -------------------------------------------------------------------------
% % Genus 0
% % -------------------------------------------------------------------------
% lgds = {};
% yy = stats(2:end, MIQ_COL);
% inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
% xx1 = xx(inds);
% yy1 = yy(inds);
% figure; hold on
% plot(xx1, cellfun(get_y_field, yy1), 'X')
% lgds(end+1) = stats(1, MIQ_COL);
% 
% yy = stats(2:end, IOQ_CONN_G0_COL);
% inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
% xx1 = xx(inds);
% yy1 = yy(inds);
% figure
% plot(xx1, cellfun(get_y_field, yy1), 'X')
% lgds(end+1) = stats(1, IOQ_CONN_G0_COL);
% 
% yy = stats(2:end, IOQ_COT_G0_COL);
% inds = cellfun(@(x) is_valid(x) && is_genus0(x), yy);
% xx1 = xx(inds);
% yy1 = yy(inds);
% figure
% plot(xx1, cellfun(get_y_field, yy1), 'X')
% lgds(end+1) = stats(1, IOQ_COT_G0_COL);
% 
% legend(lgds, 'Interpreter', 'none')
% title('Timing (genus 0)')
% hold off
% filename = fullfile(out_folder_plots, 'elapsed_total_genus0');
% print(gcf, filename, '-dpng', '-r300')
% saveas(gcf, filename, 'fig')
% 
% % -------------------------------------------------------------------------
% % highgenus
% % -------------------------------------------------------------------------
% lgds = {};
% yy = stats(2:end, MIQ_COL);
% inds = cellfun(@(x) is_valid(x) && x.genus>0, yy);
% xx1 = xx(inds);
% yy1 = yy(inds);
% try
%     close(2); 
% catch
% end
% figure(2); hold on
% plot(xx1, cellfun(get_y_field, yy1), 'X')
% lgds(end+1) = stats(1, MIQ_COL);
% 
% yy = stats(2:end, OPT1A_CONN_COL);
% inds = cellfun(@(x) is_valid(x) && x.genus>0, yy);
% xx2 = xx(inds);
% yy2 = yy(inds);
% figure(2); clf(1); hold on
% plot(xx2, cellfun(get_y_field, yy2), 'X')
% lgds(end+1) = stats(1, OPT1A_CONN_COL);
% 
% yy = stats(2:end, 7);
% inds = cellfun(@(x) is_valid(x) && x.genus>0, yy);
% xx2 = xx(inds);
% yy2 = yy(inds);
% figure(2); clf(1); hold on
% plot(xx2, cellfun(get_y_field, yy2), 'X')
% lgds(end+1) = stats(1, 7);
% 
% % yy = stats(2:end, 7);
% % inds = cellfun(@(x) is_valid(x) && x.genus>0, yy);
% % xx2 = xx(inds);
% % yy2 = yy(inds);
% % figure(2); clf(1); hold on
% % plot(xx2, cellfun(get_y_field, yy2), 'X')
% % lgds(end+1) = stats(1, 7);
% 
% % yy = stats(2:end, 8);
% % inds = cellfun(@(x) is_valid(x) && x.genus>0, yy);
% % xx2 = xx(inds);
% % yy2 = yy(inds);
% % figure(2); clf(1); hold on
% % plot(xx2, cellfun(get_y_field, yy2), 'X')
% % lgds(end+1) = stats(1, 8);
% 
% legend(lgds, 'Interpreter', 'none')
% title('Timing (genus > 0)')
% hold off
% filename = fullfile(out_folder_plots, 'elapsed_total_highg');
% print(gcf, filename, '-dpng', '-r300')
% saveas(gcf, filename, 'fig')
% 
% % plot
% % n_files = size(stats, 1) - 1;
% % n_exp = size(stats, 2) - 1;
% % legends = {};
% % figure
% % hold on
% % for c = 2:n_exp+1
% %     yy = stats(2:end, c);
% %     
% %     valid = cellfun(@(x) is_valid(x), yy);
% %     xxv = xx(valid);
% %     yyv = yy(valid);
% % 
% %     %plot(xxv, cellfun(get_y_field, yyv), ['--', ALG_TO_COLOR(stats{1, c})])
% %     plot(xxv, cellfun(get_y_field, yyv), '--')
% %     legends{end+1} = stats{1, c};
% % end
% % 
% % A = axis();
% % for i = 1:numel(ALG_OPTS)
% %     scatter(-1000, -1000, SYMBOLS(ALG_OPTS{i}), 'k', 'filled')
% %     legends{end+1} = ALG_OPTS{i};
% % end
% % axis(A);
% % 
% % label_points = false;
% % dx = 0.1;
% % dy = 0.1;
% % for c = 2:n_exp+1
% %     yy = stats(2:end, c);
% %     valid = cellfun(@(x) is_valid(x), yy);
% %     for i = 1:numel(xx)
% %         if ~valid(i)
% %             continue
% %         end
% %         
% %         scatter(xx(i), ...
% %             get_y_field(yy{i}), ...
% %             SYMBOLS(yy{i}.alg_used), ...
% %             'k', ...
% %             'filled')
% %         
% %         if label_points
% %             text(xx(i) + dx, get_y_field(yy{i}) + dy, stats{i+1, 1}, 'Interpreter', 'none');
% %         end
% %     end
% % end
% % 
% % hold off
% % title('Energy')
% % legend({legends{:}}, 'Interpreter', 'none', 'Location', 'northwest')





