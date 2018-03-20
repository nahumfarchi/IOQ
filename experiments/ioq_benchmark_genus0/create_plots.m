experiments_setup;

% x axis (number of faces)
miq_stats = stats(2:end, 4);
xx = cellfun(@(x) x.nF, miq_stats);
[xx, inds] = sort(xx);
stats_sorted = stats(2:end, 1:end);
stats_sorted = stats_sorted(inds, :);
stats(2:end, 1:end) = stats_sorted;

% y axis function
get_y_field = @(x) x.miq_energy;

is_valid = @(x) ~isempty(x) && ...
                 (~isfield(x, 'failed') || ...
                 ~x.failed);

% plot
n_files = size(stats, 1) - 1;
n_exp = size(stats, 2) - 1;
legends = {};
figure
hold on
for c = 2:n_exp+1
    yy = stats(2:end, c);
    
    valid = cellfun(@(x) is_valid(x), yy);
    xxv = xx(valid);
    yyv = yy(valid);

    plot(xxv, cellfun(get_y_field, yyv), ['--', ALG_TO_COLOR(stats{1, c})])
    legends{end+1} = stats{1, c};
end

A = axis();
for i = 1:numel(ALG_OPTS)
    scatter(-1000, -1000, SYMBOLS(ALG_OPTS{i}), 'k', 'filled')
    legends{end+1} = ALG_OPTS{i};
end
axis(A);

label_points = true;
dx = 0.1;
dy = 0.1;
for c = 2:n_exp+1
    yy = stats(2:end, c);
    valid = cellfun(@(x) is_valid(x), yy);
    for i = 1:numel(xx)
        if ~valid(i)
            continue
        end
        
        scatter(xx(i), ...
            get_y_field(yy{i}), ...
            SYMBOLS(yy{i}.alg_used), ...
            ALG_TO_COLOR(yy{i}.exp.name), ...
            'filled')
        
        if label_points
            text(xx(i) + dx, get_y_field(yy{i}) + dy, stats{i+1, 1}, 'Interpreter', 'none');
        end
    end
end

hold off
title('Energy')
legend({legends{:}}, 'Interpreter', 'none', 'Location', 'northwest')
