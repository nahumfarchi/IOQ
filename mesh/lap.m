function [lapop, W] = lap(m, E)
% LAP Return the cotangent weights laplacian of the given mesh.
%   m - the mesh
%   E - gradient matrix
% Example:
%   m = Mesh();
%   m.loadTM(path_to_off_file);
%   [gradop, E] = grad(m);
%   lapop = lap(m, E);
% or
%   lapop = lap(m);
    if nargin < 2
        [~, E] = grad(m);
    end
    W = 0.25*E'*m.Gf_inv*E;
    lapop = m.Gv_inv*W;
    %lapop = 0.25*m.Gv_inv*E'*m.Gf_inv*E;
end

