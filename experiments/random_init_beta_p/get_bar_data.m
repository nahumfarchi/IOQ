function [errY, y] = get_bar_data(data)
    errY = [];
    y = [];

    % miq energy
    X = cellfun(@(x) x.Emiq_hist(end), data);
    y(end+1) = mean(X);
    errY(end+1) = std(X);

    % E1
    X = cellfun(@(x) x.E_hist1(end), data);
    y(end+1) = mean(X);
    errY(end+1) = std(X);

    % E2
    X = cellfun(@(x) x.E_hist2(end), data);
    y(end+1) = mean(X);
    errY(end+1) = std(X);

    % E
    X = cellfun(@(x) x.E_hist(end), data);
    y(end+1) = mean(X);
    errY(end+1) = std(X);

    % m
    X = cellfun(@(x) x.m_hist(end), data);
    y(end+1) = mean(X);
    errY(end+1) = std(X);