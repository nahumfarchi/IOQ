function [N] = row_norm_square(X)
    N = sum(X.^2, 2);
end

