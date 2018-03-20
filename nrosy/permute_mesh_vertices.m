function [pm, ind] = permute_mesh_vertices(mesh)
%PERMUTE_MESH_VERTICES


ind = randperm(mesh.nV);
V = zeros(mesh.nV, 3);
V(ind, :) = mesh.V;
F = zeros(mesh.nF, 3);
for fid = 1:mesh.nF
    F(fid, 1) = ind(mesh.F(fid, 1));
    F(fid, 2) = ind(mesh.F(fid, 2));
    F(fid, 3) = ind(mesh.F(fid, 3));
end

pm = Mesh(V, F);

end

