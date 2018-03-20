%% quick clean up
% filepaths = get_filepaths(out_folder_ffields, '.ffield');
% n_files = numel(filepaths);
% progressbar
% for i = 1:n_files
%     fp = filepaths{i};
%     m = Mesh(fp);
%     if contains(fp, 'genus0') && m.genus > 0
%         delete(fp);
%     end
%     progressbar(i / n_files)
% end

%% rename files
% filepaths = get_filepaths(out_folder_ffields, '.ffield');
% n_files = numel(filepaths);
% progressbar
% for i = 1:n_files
%     fp = filepaths{i};
%     if contains(fp, 'JL')
%         delete(fp)
%     end
% %     if contains(fp, 'JL_IOQ_conn_opt1a')
% %         [p, n, e] = fileparts(fp);
% %         j = strfind(n, '_');
% %         j = j(1) - 1;
% %         n = [n(1:j), '_JL_IOQ_eps0.5.ffield'];
% %         newfp = fullfile(p, n);
% %         copyfile(fp, newfp);
% %         delete(fp);
% %     end
% %     progressbar(i / n_files);
% end

%% ========================================================================
%  Energy per edge plots
%  x axis - edge energy
%  y axis - % of edges with energy < x
%  ========================================================================
ioq_benchmark_setup;
plotting_defaults(30, 30);
print_header(experiment_name);

SAVE = true;
LOAD = true;

stats_fp = fullfile(base_folder, 'stats.mat');
if exist(stats_fp, 'file')
    stats = load(stats_fp);
    stats = stats.stats;
else
    error('stats file does not exist')
end

filepaths = get_filepaths(out_folder_ffields, '.ffield');
exp_names = stats(1, 2:end);
n_exps = numel(exp_names);
n_files = numel(filepaths);
name_to_col = containers.Map(exp_names, (1:length(exp_names))+1);
name_to_filepaths = containers.Map();
for i = 1:n_exps
    name = exp_names{i};
    %inds = cellfun(@(x) strfind(x, name), filepaths, 'UniformOutput', false);
    %inds = cellfun(@(x) ~isempty(x), inds);
    inds = cellfun(@(x) contains(x, name), filepaths);
    name_to_filepaths(name) = filepaths(inds);
end


% -------------------------------------------------------------------------
% Calculate per edge energy
% -------------------------------------------------------------------------
if ~LOAD
    name_to_edge_miq_energy = containers.Map(exp_names, cell(n_exps, 1));
    %name_to_edge_tcods_energy = containers.Map(exp_names, cell(n_exps, 1));
    %name_to_tcods_energy = containers.Map(exp_names, cell(n_exps, 1));
    name_to_miq_energy = containers.Map(exp_names, cell(n_exps, 1));
    progressbar(0, 0)
    for i = 1:n_exps
        name = exp_names{i};
        print_header(name);
        exp_filepaths = name_to_filepaths(name);
        progressbar([], 0)
        for j = 1:numel(exp_filepaths)
            fp = exp_filepaths{j};
            disp(fp)
            m = Mesh(fp);
            E = per_edge_energy(m);

            name_to_edge_miq_energy(name) = [name_to_edge_miq_energy(name); E];
            %if ~isfield(m, 'connection')
            %    disp('Creating connection...')
            %    m = TCODS(m, 'sing', [m.vert_sing; m.gen_sing]);
            %end
            %x = m.connection;
            %name_to_edge_tcods_energy(name) = ...
            %    [name_to_edge_tcods_energy(name); ...
            %     x.^2];
            % this is definitely confusing but the miq_energy field is
            % actually norm(x)^2, where x is the connection. Our assumption
            % was that these are the same.
            %name_to_tcods_energy(name) = [name_to_tcods_energy(name); norm(x)^2];
            name_to_miq_energy(name) = [name_to_miq_energy(name); sum(E)];
            
            progressbar([], j / numel(exp_filepaths))
        end
        progressbar(i / n_exps)
    end
else
    filename = fullfile(out_folder_plots, 'workspace');
    load(filename);
end

% Test that the sum of all edge energies is the same as the sum of all 
% MIQ energies
% for i = 1:n_exps
%     name = exp_names{i};
%     E_edge = name_to_edge_miq_energy(name);
%     E_mesh = name_to_tcods_energy(name);
%     assert(abs(sum(E_edge) - sum(E_mesh)) < 1e-7)
% end

if ~LOAD && SAVE
    filename = fullfile(out_folder_plots, 'workspace');
    save(filename, 'name_to_edge_miq_energy', 'name_to_miq_energy');
end

%% Plot percentage of edges with energy < x
plotting_defaults(10, 10, 1);
set(0,'DefaultFigureColormap',cbrewer('seq','YlOrRd',8));
%bins = 1000;
%bins = [0, 0.001, 0.05];
bins = 0:0.001:0.05;
lgd = {};
figure; hold on
%groups = {{'MIQ'}, {'IOQ_conn_genus0', 'IOQ_conn_opt1a'}, ...
%    {'IOQ_cot_genus0', 'IOQ_cot_opt1a'}, ...
%    {'IOQ_conn_opt2z'}, ...
%    {'IOT_cot_opt2z'}, ...
%    {'JL_IOQ_conn_opt1a_eps0.5'}};
groups = containers.Map(...
    {'MIQ', 'IOQ_conn_1a', 'JL_IOQ_eps0.5', 'GO'}, ...
    {{'MIQ'}, {'IOQ_conn_1a'}, ...
    {'JL_IOQ_eps0.5'}, ...
    {'GO'}});
for key = keys(groups)
    fprintf('key = %s\n', key{1});
    E_edges_group = [];
    for name = groups(key{1})
        fprintf('\tname = %s\n', name{1});
        E_edge = name_to_edge_miq_energy(name{1});
        fprintf('\tsum(E_edge) = %.4g\n', sum(E_edge));
        E_edges_group = [E_edges_group; E_edge];
    end
    if isempty(E_edges_group)
        continue
    end
    
    [hh, xx] = hist(E_edges_group, bins);
    yy = cumsum(hh) / length(E_edges_group) * 100;
    plot(xx, yy)
    lgd{end+1} = key{1};
end
%for i = 1:n_exps
%    name = exp_names{i};
%    E_edge = name_to_edge_miq_energy(name);
%    if isempty(E_edge)
%        continue
%    end
%    
%    [hh, xx] = hist(E_edge, bins);
%    yy = cumsum(hh) / length(E_edge) * 100;
%    plot(xx, yy)
%    lgd{end+1} = name;
%end
title('Per edge energy')
xlabel('Edge energy')
ylabel('% of edges that are < x')
%xlim([0, 0.1])
%ylim([86, 100])
legend(lgd, 'Interpreter', 'none')
if SAVE
    filename = fullfile(out_folder_plots, 'per_edge_energy');
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

%% Plot percentage of meshes with energy < x
bins = 100;
lgd = {};
figure; hold on
for key = keys(groups)
    fprintf('key = %s\n', key{1});
    E_group = [];
    for name = groups(key{1})
        fprintf('\tname = %s\n', name{1});
        E = name_to_miq_energy(name{1});
        %fprintf('\tsum(E_edge) = %.4g\n', sum(E_edge));
        E_group = [E_group; E];
    end
    if isempty(E_group)
        continue
    end
    
    [hh, xx] = hist(E_group, bins);
    yy = cumsum(hh) / length(E_group) * 100;
    plot(xx, yy)
    lgd{end+1} = key{1};
end
title('Per mesh energy')
xlabel('Mesh energy')
ylabel('% of meshes that are < x')
%xlim([0, 0.1])
%ylim([86, 100])
legend(lgd, 'Interpreter', 'none')
hold off
if SAVE
    filename = fullfile(out_folder_plots, 'per_mesh_energy');
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

return

%% ========================================================================
%  energy / time plots
%  x axis - mesh name
%  y axs - energy / time
%  ========================================================================
ioq_benchmark_setup;
plotting_defaults(30, 30, 1);
set(0,'DefaultAxesColorOrder',cbrewer('qual','Set1',8));
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
skip = {'MIQ', false; ...
    'IOQ_conn_genus0', false; ...
    'GO', false; ...
    'IOQ_conn_1a', true; ...
    'JL_IOQ_eps0.3', true; ...
    'JL_IOQ_eps0.5', true; ...
    'JL_IOQ_eps0.8', false};
    %'IOQ_cot_genus0', true; ...
    %'IOQ_conn_opt1a', false; ...
    %'IOQ_cot_opt1a', true; ...
    %'IOQ_conn_opt2z', true; ...
    %'IOQ_cot_opt2z', true; ...
    %'JL_IOQ_eps0.5', false; ...
    %'JL_IOQ_eps0.8', false; ...
    %'JL_IOQ_eps1', false; ...
    %'GO', false};
skip = containers.Map(skip(:, 1), skip(:, 2));

% sort by number of faces
miq_stats = stats(2:end, MIQ_COL);
xx = cellfun(@(x) x.nF, miq_stats);
[xx, inds] = sort(xx);
stats_sorted = stats(2:end, 1:end);
stats_sorted = stats_sorted(inds, :);
stats(2:end, 1:end) = stats_sorted;

is_valid = @(x) ~isempty(x) && ...
                 (~isfield(x, 'failed') || ...
                 ~x.failed);
             
is_genus0 = @(x) x.genus == 0;

SAVE = true;
X_FONT_SIZE = 14;

%% Plot energy
figure; hold on
lgd = {};
for i = 1:numel(exp_names)
    name = exp_names{i};
    if skip(name)
        continue
    end
    yy = stats(2:end, name_to_col(name));
    inds = cellfun(@(x) is_valid(x), yy);
    xx = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
    yy = cellfun(@(x) get_energy(x), yy);
    yy(~inds) = nan;
    if any(~isnan(yy))
        disp(['plotting ', name])
        plot(1:length(xx), yy, 'X')
        xticks(1:length(xx))
        xticklabels(xx)
        xtickangle(90)
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', X_FONT_SIZE)
        lgd{end+1} = name;
    end
end
legend(lgd, 'Location', 'northwest', 'Interpreter', 'none')
title('MIQ Energy')
hold off
if SAVE
    filename = fullfile(out_folder_plots, 'energy_xnames');
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

%% Plot timing
figure; hold on
lgd = {};
for i = 1:numel(exp_names)
    name = exp_names{i};
    if skip(name)
        continue
    end
    yy = stats(2:end, name_to_col(name));
    inds = cellfun(@(x) is_valid(x), yy);
    xx = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
    yy = cellfun(@(x) get_time(x), yy);
    yy(~inds) = nan;
    if any(~isnan(yy))
        plot(1:length(xx), yy, 'X')
        xticks(1:length(xx))
        xticklabels(xx)
        xtickangle(90)
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', X_FONT_SIZE)
        lgd{end+1} = name;
    end
end
legend(lgd, 'Location', 'northwest', 'Interpreter', 'none')
title('Timing')
hold off
if SAVE
    filename = fullfile(out_folder_plots, 'timing_xnames');
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

%% Plot singularities
figure; hold on
lgd = {};
for i = 1:numel(exp_names)
    name = exp_names{i};
    if skip(name)
        continue
    end
    yy = stats(2:end, name_to_col(name));
    inds = cellfun(@(x) is_valid(x), yy);
    xx = cellfun(@(x) strrep(x, '_', '\_'), stats(2:end, 1), 'UniformOutput', false);
    yy = cellfun(@(x) get_sing(x), yy);
    yy(~inds) = nan;
    if any(~isnan(yy))
        plot(1:length(xx), yy, '--X')
        xticks(1:length(xx))
        xticklabels(xx)
        xtickangle(90)
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', X_FONT_SIZE)
        lgd{end+1} = name;
    end
end
legend(lgd, 'Location', 'northwest', 'Interpreter', 'none')
title('Singularities')
hold off
if SAVE
    filename = fullfile(out_folder_plots, 'sing_xnames');
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

%%
S = stats;
S(2:end, 2:end) = cellfun(@get_sing, stats(2:end, 2:end), 'UniformOutput', false);
E = stats;
E(2:end, 2:end) = cellfun(@get_energy, stats(2:end, 2:end), 'UniformOutput', false);