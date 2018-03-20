FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = false;
SAVE = false;
LOAD = false;
OUT_FOLDER = 'results';
mkdir(OUT_FOLDER);

fp = '../../../data/torus_s0.off';
%fp = '../../../data/torus_fat_r2.off';
%fp = '../../../data/elephant_r.off';
%fp = '../../../data/basic/cat.off';
%fp = '../../../data/genus1_small/fertility.off';

[~, meshname, ~] = fileparts(fp);
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;
alpha_g = get_gaussian_curvature(m);
beta_g = wrapToPi(generator_angle_defects(m));

% run IOQ
rng(SEED)
[alpha, beta, elapsed_ioq] = IOQ_highgenus_gpu(...
    V, F, ...
    'UseGPU', USE_GPU, ...
    'Iterations', 2000);

k = [alpha; beta];
res_ioq = TCODS(m, ...
    'k', k, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);

d0 = get_exterior_derivatives(m);
L = d0'*d0;
H = -m.H;
%H = [H(:, 1), -H(:, 2)];
%H = [-H(:, 1), H(:, 2)];
F = factorize(L);
A2 = (H'*d0) / F;
y0 = (2/pi)*(A2*alpha_g - beta_g);
B = H - d0 * (F \ (d0' * H));
check_norm('round(A2*alpha - y0)', 'beta', 'Log', -1);

%M2 = B * inverse(H'*B);
%[x, uk] = closest_lpoint_vor(full(M2*A2), M2*y0);

%X = (pi/2) * B * inverse(H'*B);
%a = inverse(L) * ((pi/2)*alpha - alpha_g);
%pt = X*(beta_g + H'*d0*a);
%[x, uk] = closest_lpoint_vor(full((pi/2)*X), pt);
%b1 = inverse(H'*B)*((pi/2)*beta-beta_g-H'*d0*a);
%b2 = inverse(H'*B)*((pi/2)*uk-beta_g-H'*d0*a);
%x = d0*a + B*b2;
%%

a1 = inverse(L) * ( (pi/2)*alpha - alpha_g );
b1 = inverse(H'*B)*((pi/2)*beta-beta_g-H'*d0*a1);
E1 = norm(d0*a1)^2 + norm(B*b1)^2;

X = B * inverse(H'*B) * (pi/2);
pt = B * inverse(H'*B) * ( beta_g + H'*d0*a1 );
X2 = X;
%X2 = [X(:,1), X(:,2)];
%X2'*X2
%[~, uk] = closest_lpoint_vor(full(X2), pt);

%a2 = inverse(L) * ( (pi/2)*alpha - alpha_g );
%b2 = inverse(H'*B)*((pi/2)*uk-beta_g-H'*d0*a1);
%E2 = norm(d0*a2)^2 + norm(B*b2)^2;

%fprintf('|d0 a1| = %.2f, |B b1| = %.2f\n', norm(d0*a1)^2, norm(B*b1)^2);
%fprintf('|d0 a2| = %.2f, |B b2| = %.2f\n', norm(d0*a2)^2, norm(B*b2)^2);
%fprintf('E1, E2 = %.2f, %.2f\n', E1, E2);

%%
%   m    = 4;                        % 4x4 MIMO system
%   xMax = 2;                        % 16-QAM modulation
%   cplx = 1;
%   y    = randn(m,1)+i(randn(m,1);  % generate random complex target vector
%   H    = randn(m,m)+i*randn(m,m);  % generate random complex channel 
%   [xHat,wHat,nv,nf,nl] = SEA_det(m,xMax,y,H,cplx);
y = pt;
m = size(X2, 2);
xMax = 20;
cplx = false;
[xHat,wHat,nv,nf,nl] = SEA_det(m,xMax,y,X2,cplx);

a3 = inverse(L) * ( (pi/2)*alpha - alpha_g );
b3 = inverse(H'*B)*((pi/2)*xHat-beta_g-H'*d0*a1);
E3 = norm(d0*a3)^2 + norm(B*b3)^2;

fprintf('%.2f, %.2f\n', E1, E3);

return

%%
k2 = [alpha; uk];
res_ioq2 = TCODS(m, ...
    'k', k2, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);
fprintf('%.2f, %.2f\n', res_ioq.miq_energy, res_ioq2.miq_energy);
fprintf('Rou |Bx-y0| = %.2f\n', norm(X*beta-pt)^2);
fprintf('Vor |Bx-y0| = %.2f\n', norm(X*uk-pt)^2);

%assert( abs( (res_ioq.miq_energy - res_ioq2.miq_energy) - ...
%             (norm(B*b2)^2 - norm(B*b1)^2) ) ...
%             < 1e-10)

%%

figure
title_ioq = {sprintf('IOQ round, E = %.2f', res_ioq.miq_energy), ...
    sprintf('|Bx-y0| = %.2f', norm(X*beta-pt));};
title_ioq2 = {sprintf('IOQ vor, E = %.2f', res_ioq2.miq_energy), ...
    sprintf('|Bx-y0| = %.2f', norm(X*uk-pt));};
subplot(121); res_ioq.draw; title(title_ioq)
subplot(122); res_ioq.draw; title(title_ioq2)
