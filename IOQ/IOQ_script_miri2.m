addpath(genpath('../data/'));
clear all
close all

run_miq = 1;

rng(112); % Fix the seed
%%%%%
FACE0 = 1;        % starting face for tcods
THETA0 = nan;     % (use global constraint instead)
GVEC = [1, 0, 0]; % constraint vector (in global coordinates)
DEGREE = 4;       % works only for degree 4 atm

disp('Loading mesh...')
name = 'cup';
mesh = Mesh([name '.off']);
V = mesh.V;
F = mesh.F;
% or: 
%mesh = Mesh(V, F);

nv = mesh.nV; nf = mesh.nF; ne = mesh.nE;


%%%%%
% disp('Running IOQ...')
[k,Lp,E_hist] = IOQ_cpu2(V, F, 'Laplacian', 'conn2','Plot',true,'Iterations',1000);
S = find(k); kk = k(S);
mesh.vert_sing = [S,kk];
figure;mesh.draw;

% Add non-contractible cycle holonomies
g = zeros((1-round(sum(k))/8)*2,1);
k = [k;g];

%%%%%
disp('Creating direction field...')
[A, K, d0, d1, H] = tcods_gsystem(V, F);

%figure(); mesh.draw(); hold on; mesh.plotH(); hold off

% P = full(H'*H);
% [uu,ss,vv] = svd(P);

n = size(A, 2); 
b = K - (pi / 2) * k;
x = lsqlin(speye(n, n), zeros(n, 1), [], [], A, -b, -inf(n, 1), inf(n, 1));

assert(norm(A*x - (-b)) < 1e-10)
assert(norm(d1*x) < 1e-10)
%assert(abs(rank(full([A; d1])) - min(size([A; d1]))) < 1e-10)

norm(x)^2
return

u = d0\x; u = u - mean(u);
fprintf('u range: %g\n', max(u) - min(u));

b1 = b(1:nv); b2 = b(nv+1:end);

% Basis for harmonic forms
B = H - d0*(d0\H);
assert(norm(d0'*B,'fro') < 1e-10)
assert(norm(d1*B,'fro') < 1e-10)


% Looking for x = d0*u1 + B*u2
% s.t. [d0';H']x + b = 0
% (1) d0'*x + b1 = 0 and (2) H'x + b2 = 0
% For (1): 
% d0'*(d0*u1 + B*u2) + b1 = 0
% d0'*d0*u1 + d0'*B*u2 + b1 = 0
% h is harmonic, so d0'*B = 0
% d0'*d0*u1 + b1 = 0 ==> u1 = (d0'*d0)\(-b1)
% For (2):
% H'*d0*u + H'*B*u2 + b2 = 0
% H'*d0*[(d0'*d0)\(-b1)] + H'*B*u2 + b2 = 0
% u2 = (H'*B)\(-b2 - H'*d0*[(d0'*d0)\(-b1)])
% u2 = (H'*B)\(-b2) - H'\(H'*d0*[(d0'*d0)\(-b1)])
u1 = (d0'*d0)\(-b1); 
u2 = (H'*B)\(-b2) - (H'*B)\(H'*d0*[(d0'*d0)\(-b1)]);
norm(x - (d0*u1 + B*u2))

% According to: https://en.wikipedia.org/wiki/Block_matrix_pseudoinverse
% x = (d0' P_H)\b1 + (H' P_d0)\b2
P_H = speye(ne) - H*inv(H'*H)*H';
P_d0 = speye(ne) - d0*pinv(full(d0));
%norm(x - (pinv(full(d0' * P_H))*(-b1) + pinv(full(H'*P_d0))*(-b2)))
%norm(x - (pinv(full(d0' * P_H))*(-b1) + pinv(full(B'))*(-b2)))
norm(x - (P_H*d0*pinv(full(d0'*P_H*d0))*(-b1) + pinv(full(B'))*(-b2)))
norm(b1'*(pinv(full(d0' * P_H))'*pinv(full(d0' * P_H)))*b1 - b1'*pinv(full(d0'*P_H*d0))*b1)

fprintf('%g %g %g\n',norm(x)^2, b1'*pinv(full(d0'*P_H*d0))*b1, b2'*pinv(full(B'))'*pinv(full(B'))*b2)

%M1 = d0' - d0'*H*inv(H'*H)*H';
%norm(pinv(full(M1))-P_H*d0*pinv(full(d0'*P_H*d0)),'fro')

% Min ||x||^2 = ||d0*u1||^2 + ||B*u2||^2
% = ||d0*pinv(d0'*d0)*b1||^2 + ||B*inv(H'*B)*(b2 - H'*d0*pinv(d0'*d0)*b1)||^2
%%% gr = round((K(nv+1:end)-H'*d0*u1)/(pi/2));
%gr = round((K(nv+1:end) + H'*d0*Lp*(-b1))/(pi/2));
% kr = [k(1:nv);gr];
% br = K - (pi / 2) * kr;
% b1r = br(1:nv); b2r = br(nv+1:end);
% b2r = -H'*d0*u1;
% b2r = K(nv+1:end) - (pi / 2) * gr;
%xr = lsqlin(speye(n, n), zeros(n, 1), [], [], A, -br, -inf(n, 1), inf(n, 1));

%u2r = (H'*B)\(-b2r - H'*d0*u1);
b_h = K(nv+1:end) - H'*d0*Lp*K(1:nv) + ...
    H'*d0*Lp*(pi / 2) * k(1:nv);
%gr = round(b_h/(pi/2));
u2r = inv(H'*B)*((pi / 2) * round(b_h/(pi/2)) - b_h);
e_round = norm(B*u2r)^2;

% CVP for gr
% argmin_{g_r \in Z} ||B*inv(H'*B)*(gr - b_h/(pi/2))||^2
% Bruteforce for small genus
hh = length(K)-nv;
bs = 4;
e = zeros(bs^hh,1);
for i=0:bs^hh-1
    gr = de2bi(i,hh,bs)';
    e(i+1) = norm(B*inv(H'*B)*((pi/2)*gr - b_h))^2;
end
[min(e), e_round]
norm(d0*u1)^2 + [min(e), e_round]
xr = d0*u1 + B*u2r;

fprintf('E IOQ : %g %g\n', norm(x)^2, norm(xr)^2);

% u1r = d0\xr; norm(u1 - mean(u1) - (u1r - mean(u1r)))
% 
% b1r = br(1:nv); b2r = br(nv+1:end);
% u1r = (d0'*d0)\(-b1r); 
% u2r = (H'*B)\(-b2r - H'*d0*((d0'*d0)\(-b1r)));
% norm(xr - (d0*u1r + B*u2r))
% norm(xr)^2 - (norm(d0*u1r)^2 + norm(B*u2r)^2)
% norm((-b2r - H'*d0*((d0'*d0)\(-b1r))))
% 


[ffield] = connection_to_nrosy(...
                mesh, ...
                x, ...
                FACE0, ...
                THETA0, ...
                DEGREE, ...
                'gConstraintVec', GVEC);


%%%%%
disp('Plotting...')
mesh.ffield_vectors = ffield;
mesh.degree = DEGREE;
S = find(k(1:mesh.nV)); kk = k(S);
mesh.vert_sing = [S,kk];
figure()
mesh.draw()

disp('Running MIQ param...')
R = mesh.ffield_vectors(1:nf, :);
[uv, fuv] = MIQ_param_mex_bin(mesh.V, mesh.F, R, 75);

disp('Saving...')
%options.nm_file = ['ext/rosy2gridparam' 'cross.png'];
options.nm_file = './cross.png';
options.face_texcorrd = fuv;
options.object_texture = uv;

write_obj_for_quad('../data/', ...
    [name '_uv.obj'], ...
    mesh.V, ...
    mesh.F, ...
    options);

disp('Saving quad...')
system(['..\..\..\..\..\software\rosy2gridparam_mex\external\QEX_bin', ' ', '../data/' name '_uv.obj', ' ', '../data/' name '_quad.obj']);

% MIQ
if run_miq
    [theta, p, k_miq, R, local_frames, frame_diffs, elapsed, pFixed] = NRoSy_mex_bin(['../data/' name '.off' ], [1], [1,0,0], 4);
    fprintf('E MIQ : %g\n', E_MIQ(mesh, theta, frame_diffs, p, 4));
    S_miq = find(k_miq(1:mesh.nV)); kk = -k_miq(S_miq);
    mesh.vert_sing = [S_miq,kk];
    figure()
    mesh.draw()
end