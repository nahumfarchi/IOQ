function [edges v2e] = compute_edges(face)

% compute_edges - from a list of faces, compute the list of (unique) edges.
%
%   [edge v2e] = compute_edges(face);
%   v2e returns a sparse matrix such that v2e(a, b) gives the
%   index i in edge such that edge(1, i) = a; edge(2, i) = b.
%
%   works for triangular and tets meshes.
%
%   Copyright (c) 2004 Gabriel Peyr�

if isempty(face)
    edges=[];
    return;
end

[tmp,face] = check_face_vertex([],face);

if size(face,1)~=3 && size(face,1)~=4
    error('Problem, works for triangles and tets only.');
end

d = size(face,1);

edges = [];
for i=1:d
    sel = [i, mod(i,d)+1];
    edges = [edges, face(sel,:)];
end

% sort pair of vertex
I = find(edges(1,:)>edges(2,:));
J = find(edges(1,:)<=edges(2,:));
edges = [edges(end:-1:1,I), edges(:,J)];

% unique id
m = max(edges(:))+100;
id = edges(1,:) + m*edges(2,:);

[tmp,I] = unique(id);
edges = edges(:,I);

% added by dilip to compute v2e
ne = length(edges);
v2e = sparse([edges(1, :)' edges(2, :)'], [edges(2, :)' edges(1, :)'], [1:ne 1:ne]');


