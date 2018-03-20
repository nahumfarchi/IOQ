function [edgelist] = create_dual_edgelist(mesh)
%CREATE_EDGELIST Create an edgelist from the given mesh (of the dual graph).

edgelist = zeros(mesh.niE, 2);
EF = mesh.EFAdj;

row = 1;
for eid = 1:mesh.nE
    di = EF(eid, 1);
    dj = EF(eid, 2);
    if di < 1 || dj < 1
        warning('Boundary edges are not well tested yet.')
        continue;
    end
    
%     assert(di < dj);
    
    edgelist(row, :) = [di, dj];
    row = row + 1;
end

end

