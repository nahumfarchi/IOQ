function t = get_time(x)
    if isfield(x, 'elapsed_total')
        t = x.elapsed_total;
    else
        t = nan;
    end