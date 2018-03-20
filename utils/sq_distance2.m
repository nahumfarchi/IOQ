function D = sq_distance2(X, Y)
% Pairwise distance between columns
D = bsxfun(@plus,dot(X,X,1)',dot(Y,Y,1))-2*(X'*Y);