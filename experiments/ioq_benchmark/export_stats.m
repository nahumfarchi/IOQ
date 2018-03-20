ioq_benchmark_setup;
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
    'IOQ_cot_genus0', true; ...
    'IOQ_conn_opt1a', false; ...
    'IOQ_cot_opt1a', true; ...
    'IOQ_conn_opt2z', true; ...
    'IOQ_cot_opt2z', true; ...
    'JL_IOQ_eps0.5', false};
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

n_files = size(stats, 1) - 1;
j1 = name_to_col('IOQ_conn_genus0');
j2 = name_to_col('IOQ_conn_opt1a');
for i = 2:n_files+1  
    st1 = stats{i, j1};
    st2 = stats{i, j2};
    if is_valid(st1) && ~is_genus0(st1)
        stats(i, j1) = stats(i, j2);
    elseif is_valid(st2) && ~is_genus0(st2)
        stats(i, j1) = stats(i, j2);
    end
end

stats_clean = cell(n_files+1, 4);
stats_clean(:, 1:3) = stats(:, 1:3);
stats_clean(:, 4) = stats(: ,end);
stats_clean{1, 3} = 'IOQ_conn';
filename = fullfile(base_folder, 'stats_clean');
save(filename, 'stats_clean');
