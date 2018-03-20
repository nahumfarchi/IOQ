function [E_edges, Emiq] = per_edge_energy(mesh)
% [E_edges, Emiq] = per_edge_energy(mesh)

tol = 1e-7;

ne = mesh.nE;
E_edges = zeros(ne, 1);
p = zeros(ne, 1);
t = mesh.ffield_angles;
k = mesh.frame_diffs;

EF = mesh.EFAdj;
for eid = 1:ne
    i = EF(eid, 1);
    j = EF(eid, 2);
    ti = t(i);
    tj = t(j);
    kij = k(eid);
    p(eid) = round((tj - ti - kij) / (pi/2));
    E_edges(eid) = (ti + kij + (pi/2)*p(eid) - tj)^2;
end

if nargout > 1
    Emiq = sum(E_edges);
end

assert(abs(mesh.miq_energy - sum(E_edges)) < tol)
if ~isempty(mesh.connection)
    x = mesh.connection;
    assert(abs(norm(x)^2 - sum(E_edges)) < tol)
end

end

