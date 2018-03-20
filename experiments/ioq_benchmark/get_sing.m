function s = get_sing(x)
    if isfield(x, 'miq_energy')
        s = x.n_singularities;
    else
        s = nan;
    end