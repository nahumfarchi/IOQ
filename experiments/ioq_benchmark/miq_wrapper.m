function res = miq_wrapper(fp, varargin)
    [res.m, res.elapsed_total] = nrosy_mex(fp, varargin{:});
    res.alg_used = 'miq';
end

