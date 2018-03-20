function D = sq_distance(X)
% Pairwise distance between columns
tmp = dot(X, X);
D = bsxfun(@plus, tmp', tmp - 2*(X'*X));
%D = bsxfun(@plus, dot(X,X,1)',dot(X,X,1))-2*(X'*X);