function [N] = row_norm(X)
    N = sqrt(sum(X.^2, 2));
end

