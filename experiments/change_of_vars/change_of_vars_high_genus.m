%% =======================================================================
%  Change of vars on high genus.
%  =======================================================================
clear all
close all

FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = false;

fp = '../../../data/bunny_simple.off';
m = Mesh(fp);

%% Set alpha, beta using IOQ
V = m.V; F = m.F; ne = m.nE; nv = m.nV; nf = m.nF; ng = m.genus;
% Use IOQ
rng(SEED);
[alpha, beta] = IOQ_highgenus_gpu(V, F, 'UseGPU', USE_GPU, 'highg_method', 'option1a');

%% Geometry stuff
H = -m.H;
[d0, d1] = get_exterior_derivatives(m);
[frames, r] = create_local_frames(m);
alpha_g = get_gaussian_curvature(m);
beta_g = wrapToPi(generator_angle_defects(m));

%% Setup the system
alpha_h = alpha_g - (pi/2) * alpha;
beta_h = beta_g - (pi/2) * beta;

A = [d0'; H'];
b = (2 / pi) * ([alpha_h; beta_h] - A*r);

%% TUM mystery
% [ii,jj] = find(H); ilocs = setdiff(1:size(H,1),ii); assert(norm(H(ilocs,:),'fro')==0);
% %tic; assert(is_TUM(H(ii,:))); toc
locs = find(sum(abs(A))>=3);
Ac = A(:,locs);
% % One row from d0 and one row from H
% for j=1:2*ng
%     tic
%     for i=1:nv
%         assert(is_TUM(Ac([i,nv+j],:)))
%     end
%     toc
% end

%% Solve for all periods
p = A \ b;
% Check solution was found
assert(norm(A*p - b,'fro') < 1e-10)
% Check that p is integral
assert(norm(round(p) - p) < 1e-10)

%% Set constrained periods to 0
FF = sparse([1:nf,1:nf,1:nf]',[m.FFAdj(:,1);m.FFAdj(:,2);m.FFAdj(:,3)],ones(nf*3,1),nf,nf); assert(norm(FF-FF','fro')==0);
[it,jt,vt] = mst(FF);
EE = sparse([m.EFAdj(:,1);m.EFAdj(:,2)], [m.EFAdj(:,2);m.EFAdj(:,1)], [1:ne, -[1:ne]]', nf, nf);
p_zero_locs = abs(EE(sub2ind(size(EE),it,jt))); assert(length(p_zero_locs) == nf-1);

p_free = setdiff([1:ne],p_zero_locs);
% Pz = speye(ne,ne); Pz = Pz(p_zero_locs, :);
% bz = sparse(length(p_zero_locs),1);
% A2 = [A; Pz]; b2 = [b; bz];
A2 = A(2:end,p_free);
b2 = b(2:end);

p2 = A2 \ b2;
% Check solution was found
assert(norm(A2*p2 - b2,'fro') < 1e-10);
% Check that p is integral
assert(norm(round(p2) - p2) < 1e-10)

det(A2)

%%%%%%%% d0^T has co-rank 1 --> can remove one row as the others are
%%%%%%%% linearly dependent on it. After adding H^T and removing columns 
%%%%%%%% with fixed p values, the resulting matrix is square. Need to
%%%%%%%% show its determinant is +1 or -1 to show solution is integer and
%%%%%%%% unique.

return

%% Use random alpha, beta
for kk = 1:100
    alpha_g0 = (2/pi) * alpha_g;
    c = round(abs(sum(alpha_g0))); % = 4 xi
    n_pos_sing = (c + round(sum(alpha_g0))) / 2;
    n_neg_sing = (c - round(sum(alpha_g0))) / 2;
    inds_pos = randperm(nv, n_pos_sing);
    inds_neg = randperm(nv, n_neg_sing);
    alpha = zeros(nv, 1);
    alpha(inds_pos) = 1;
    alpha(inds_neg) = -1;
    beta = floor(4*rand(2*ng, 1)) - 2;

    %% Setup system
    alpha_h = alpha_g - (pi/2) * alpha;
    beta_h = beta_g - (pi/2) * beta;

    A = [d0'; H'];
    b = (2 / pi) * ([alpha_h; beta_h] - A*r);
    b = round(rand(size(b))); b(1) = -sum(b(2:nv));
    
    %% Solve for periods
%     p = A \ b;
%     % Check solution was found
%     assert(norm(A*p - b,'fro') < 1e-10)
%     % Check that p is integral
%     assert(norm(round(p) - p) < 1e-10)

    A2 = [A; Pz]; b2 = [b; bz];
    p2 = A2 \ b2;
    % Check solution was found
    assert(norm(A2*p2 - b2,'fro') < 1e-10);
    % Check that p is integral
    assert(norm(round(p2) - p2) < 1e-10)

end