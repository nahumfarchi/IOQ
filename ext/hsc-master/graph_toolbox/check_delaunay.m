function [ic cots] = check_delaunay(vertex, edges)

% check_delaynay - given edge (va, vb) with common neighbors vc and vd
% see if (va, vb) needs flipping. This is done as follows:
% find lengths of all 5 edges: (va, vb), (va, vc), (vb, vc), (va, vd), (vb, vd)
% then use the formulas given in the paper of Fisher et. al.'s Intrinsic Delaunay 
% Triangulation
%     va
%    / | \
%  vc  |  vd
%   \  |  /
%     vb 
% if vc or vd are 0, ignore them (boundary vertex)

va = edges(1, :);
vb = edges(2, :);
vc = edges(3, :);
vd = edges(4, :);

idxc = find(vc > 0);
idxd = find(vd > 0);
nnzc = vc > 0;
nnzd = vd > 0;

vab = find_length(vertex, va, vb);
vac = zeros(size(vab));
vbc = zeros(size(vab));
vad = zeros(size(vab));
vbd = zeros(size(vab));
vac(:, idxc) = find_length(vertex, va(:, idxc), vc(:, idxc));
vbc(:, idxc) = find_length(vertex, vb(:, idxc), vc(:, idxc));
vad(:, idxd) = find_length(vertex, va(:, idxd), vd(:, idxd));
vbd(:, idxd) = find_length(vertex, vb(:, idxd), vd(:, idxd));

% angle acb
tan_acb = find_angle(vab, vac, vbc);
% angle adb
tan_adb = find_angle(vab, vad, vbd);

cot_acb = (1 - tan_acb.^2)./(2*tan_acb);
cot_adb = (1 - tan_adb.^2)./(2*tan_adb);

cots = nnzc.*cot_acb + nnzd.*cot_adb;
ic = (cots > 0);

%lens{1} = vab;
%lens{2} = vac;
%lens{3} = vbc;
%lens{4} = vad;
%lens{5} = vbd;

function lens = find_length(vertex, va, vb)
	lens = sqrt(sum((vertex(:, va) - vertex(:, vb)).^2, 1));

% given a triangle with sides of lenth a, b, and c
% return tan (alpha/2) where alpha is the angle opposite
% side of length a
function tan_alpha2 = find_angle(a, b, c)
	tan_alpha2 = sqrt(((a-b+c).*(a+b-c))./((a+b+c).*(-a+b+c)));
