function res = go_wrapper(fp, varargin)
    [res.m, res.elapsed_total] = GO(fp, varargin{:});
    res.alg_used = 'go';
end

