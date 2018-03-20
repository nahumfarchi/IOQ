function plot_data(data)
    % miq energy
    hold on
    X = cellfun(@(x) x.Emiq_hist, data, 'UniformOutput', false);
    X = padcat(X{:});
    n = size(X, 2);
    xx = 1:n;
    yy = nanmean(X, 1);
    ee = nanstd(X);
    errorbar(xx, yy, ee);

    % E1
    X = cellfun(@(x) x.E_hist1, data, 'UniformOutput', false);
    X = padcat(X{:});
    n = size(X, 2);
    xx = 1:n;
    yy = nanmean(X, 1);
    ee = nanstd(X);
    errorbar(xx, yy, ee);

    % E2
    X = cellfun(@(x) x.E_hist2, data, 'UniformOutput', false);
    X = padcat(X{:});
    n = size(X, 2);
    xx = 1:n;
    yy = nanmean(X, 1);
    ee = nanstd(X);
    errorbar(xx, yy, ee);

    % E
    X = cellfun(@(x) x.E_hist, data, 'UniformOutput', false);
    X = padcat(X{:});
    n = size(X, 2);
    xx = 1:n;
    yy = nanmean(X, 1);
    ee = nanstd(X);
    errorbar(xx, yy, ee);

    % m
    X = cellfun(@(x) x.m_hist, data, 'UniformOutput', false);
    X = padcat(X{:});
    n = size(X, 2);
    xx = 1:n;
    yy = nanmean(X, 1);
    ee = nanstd(X);
    errorbar(xx, yy, ee);

    % gridsearch ticks
    X = cellfun(@(x) x.gridsearch_ticks, data, 'UniformOutput', false);
    X = any(padcat(X{:}) > 0, 1);
    idx = find(X);
    xx = idx;
    tmp = ylim;
    yy = ones(1, length(xx))*tmp(1);
    %X = cellfun(@(x) find(x), X, 'UniformOutput', false);
    plot(xx, yy, 'rX')
    
    hold off
    legend('Emiq', 'E1', 'E2', 'E', 'm', 'gridsearch ticks')