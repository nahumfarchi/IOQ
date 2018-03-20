function ne = get_ne(x)
    if isfield(x, 'nE')
        ne = x.nE;
    else
        ne = nan;
    end