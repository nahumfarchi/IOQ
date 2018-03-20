%% =======================================================================
%  Change of vars on high genus.
%  MIQ(mesh) --> (theta1, p1, r1) --> (x1, alpha1, beta1) --> (theta2, p2)
%  Check that E(theta1, p1) == |x1|^2 == E(theta2, p2)

%  TCODS(alpha1, beta1) --> x2 --> (theta3, p3)
%  Check that E(theta1, p1) == |x2|^2 == E(theta3, p3)
%  =======================================================================
clear all
close all

FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = false;
LOG = -1;

%fp = '../../../data/torus_s0.off';
%fp = '../../../data/3holes_lo.off';
fp = '../../../data/3holes.off';
%fp = '../../../data/yeah_right.off';
m = Mesh(fp);
nv = m.nV; nf = m.nF; ne = m.nE;
ng2 = 2*m.genus;
H = -m.H;
[d0, d1] = get_exterior_derivatives(m);
%[frames, r] = create_local_frames(m); % use MIQ frames instead
alpha_g = get_gaussian_curvature(m);
beta_g = wrapToPi(generator_angle_defects(m));

%% MIQ --> (theta1, p1, r1)
disp('MIQ --> (theta1, p1, r1)')
[miq1, elapsed_miq] = ...
    nrosy_mex(fp, FACE0, GVEC, DEGREE);
theta1 = miq1.ffield_angles;
p1 = miq1.periods;
r1 = miq1.frame_diffs;

%% (theta1, p1, r1) --> (x1, alpha1, beta1)
disp('(theta1, p1, r1) --> (x1, alpha1, beta1)')
alpha1 = round((2/pi) * (alpha_g - d0'*r1) - d0'*p1);
beta1 = round((2/pi) * (beta_g - H'*r1) - H'*p1);
x1 = -r1 - d1'*theta1 - (pi/2)*p1;
check_norma('norm(x1)^2', 'miq1.miq_energy', 'Log', LOG);

%% (x1, alpha1, beta1) --> (theta2, p2)
disp('(x1, alpha1, beta1) --> (theta2, p2)')
alpha_h = alpha_g - (pi/2) * alpha1;
beta_h = beta_g - (pi/2) * beta1;

A = [d0'; H'];
b = (2 / pi) * ([alpha_h; beta_h] - A*r1);

FF = sparse([1:nf,1:nf,1:nf]', ...
            [m.FFAdj(:,1);m.FFAdj(:,2);m.FFAdj(:,3)], ...
            ones(nf*3,1), ...
            nf,nf); 
assert(norm(FF-FF','fro')==0);
[it,jt,vt] = mst(FF);
EE = sparse([m.EFAdj(:,1);m.EFAdj(:,2)], ...
            [m.EFAdj(:,2);m.EFAdj(:,1)], ...
            [1:ne, -[1:ne]]', ...
            nf, nf);
p_zero_locs = abs(EE(sub2ind(size(EE),it,jt))); 
assert(length(p_zero_locs) == nf-1);

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
assert(norm(round(p2) - p2) < 1e-9)
%assert(norm(A2*mod(round(p2),4)-b2) < 1e-10);

det(A2)

%%%%%%%% d0^T has co-rank 1 --> can remove one row as the others are
%%%%%%%% linearly dependent on it. After adding H^T and removing columns 
%%%%%%%% with fixed p values, the resulting matrix is square. Need to
%%%%%%%% show its determinant is +1 or -1 to show solution is integer and
%%%%%%%% unique.
%%
tmp = zeros(ne, 1);
tmp(p_free) = round(p2);
%tmp = mod(round(tmp), 4);

% check that the this works with the MIQ periods
%theta2 = d1' \ (-r1 - x1 - (pi/2)*p1);  
%check_norma('d1''*theta2', '-r1 - x1 - (pi/2)*p1', 'Log', LOG);
%check_norma('miq1.miq_energy', 'E_MIQ(m, theta2, r1, p1, DEGREE)');

theta2 = d1' \ (-r1 -x1 - (pi/2)*tmp);
check_norma('d1''*theta2', '-r1 - x1 - (pi/2)*tmp', 'Log', LOG);
check_norma('miq1.miq_energy', 'E_MIQ(m, theta2, r1, tmp, DEGREE)', 'Log', LOG);

%% Solve TCODS (alpha1, beta1)
disp('TOCDS(alpha1, beta1) --> x2')
k1 = [alpha1; beta1];
res_tc = TCODS(m, ...
    'k', k1, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);
x2 = res_tc.connection;
check_norma('miq1.miq_energy', 'res_tc.miq_energy', 'Log', LOG);
check_norma('x1', 'x2', 'Log', LOG);
% x2 comes out the same as x1 so I didn't bother to go from x2 to (theta3, p3)

% Check that that angles from integrating x are the same as theta1 (up to
% mod pi/2)
theta_intg_tc = res_tc.ffield_angles;
check_norma('mod(theta1, pi/2)', 'mod(theta_intg_tc, pi/2)', 'Log', LOG, 'Tol', 1e-9);