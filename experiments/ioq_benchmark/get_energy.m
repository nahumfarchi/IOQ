function e = get_energy(x)
    if isfield(x, 'miq_energy')
        e = x.miq_energy;
    else
        e = nan;
    end