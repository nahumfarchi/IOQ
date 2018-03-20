function e = get_nF(x)
    if isfield(x, 'nF')
        e = x.nF;
    else
        e = nan;
    end