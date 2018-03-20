function [X] = normalize_rows(X)
    X = X ./ repmat(row_norm(X), 1, 3);
end

