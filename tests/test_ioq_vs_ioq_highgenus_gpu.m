%%%%%
FACE0  = 1;         % starting face for tcods
THETA0 = nan;       % (use constraint in 3D coords instead)
GVEC   = [1, 0, 0]; % constraint vector (in global coordinates)
DEGREE = 4;         % works only for degree 4 atm
SEED   = 112;
TOL    = 1e-6;
PLOT_PROPS = {'FaceColor', 'w', 'PlotField', true};

%%
m = Mesh('../data/bunny.off');
V = m.V;
F = m.F;
% or: 
%mesh = Mesh(V, F);

nv = m.nV; nf = m.nF;

rng(SEED)
[alpha1, beta1] = IOQ_highgenus_gpu(V, F);
res1 = TCODS(m, ...
             'k', [alpha1; beta1], ...
             'f0', FACE0, ...
             'theta0', THETA0, ...
             'degree', DEGREE, ...
             'CreateFField', true, ...
             'Duplicate', true, ...
             'gConstraintVec', GVEC);

rng(SEED)
[alpha2, beta2, x2, stats2, out2] = IOQ(V, F);

check_norm('res1.miq_energy', 'norm(x2)^2', 'Log', -1, 'Tol', TOL);
check_norm('alpha1', 'alpha2', 'Log', -1);

%%
m = Mesh('../data/3holes.off');
V = m.V;
F = m.F;
% or: 
%mesh = Mesh(V, F);

nv = m.nV; nf = m.nF;

rng(SEED)
[alpha1, beta1] = IOQ_highgenus_gpu(V, F);
res1 = TCODS(m, ...
             'k', [alpha1; beta1], ...
             'f0', FACE0, ...
             'theta0', THETA0, ...
             'degree', DEGREE, ...
             'CreateFField', true, ...
             'Duplicate', true, ...
             'gConstraintVec', GVEC);

rng(SEED)
[alpha2, beta2, x2, stats2, out2] = IOQ(V, F);

check_norm('res1.miq_energy', 'norm(x2)^2', 'Log', -1, 'Tol', TOL);
check_norm('[alpha1; beta1]', '[alpha2; beta2]', 'Log', -1);

%%
m = Mesh('../data/torus_s0.off');
V = m.V;
F = m.F;
% or: 
%mesh = Mesh(V, F);

nv = m.nV; nf = m.nF;

rng(SEED)
[alpha1, beta1] = IOQ_highgenus_gpu(V, F);
res1 = TCODS(m, ...
             'k', [alpha1; beta1], ...
             'f0', FACE0, ...
             'theta0', THETA0, ...
             'degree', DEGREE, ...
             'CreateFField', true, ...
             'Duplicate', true, ...
             'gConstraintVec', GVEC);

rng(SEED)
[alpha2, beta2, x2, stats2, out2] = IOQ(V, F);

check_norm('res1.miq_energy', 'norm(x2)^2', 'Log', -1, 'Tol', TOL);
check_norm('[alpha1; beta1]', '[alpha2; beta2]', 'Log', -1);

%%
m = Mesh('../data/torus_fat_r2.off');
V = m.V;
F = m.F;
% or: 
%mesh = Mesh(V, F);

H = m.H;
beta_g = m.generator_defects;

nv = m.nV; nf = m.nF;

rng(SEED)
[alpha1, beta1] = IOQ_highgenus_gpu(V, F);
res1 = TCODS(m, ...
             'k', [alpha1; beta1], ...
             'f0', FACE0, ...
             'theta0', THETA0, ...
             'degree', DEGREE, ...
             'CreateFField', true, ...
             'Duplicate', true, ...
             'gConstraintVec', GVEC);

rng(SEED)
[alpha2, beta2, x2, stats2, out2] = IOQ(V, F);
res2 = TCODS(m, ...
             'k', [alpha2; beta2], ...
             'f0', FACE0, ...
             'theta0', THETA0, ...
             'degree', DEGREE, ...
             'CreateFField', true, ...
             'Duplicate', true, ...
             'gConstraintVec', GVEC);

check_norm('res1.miq_energy', 'norm(x2)^2', 'Log', -1, 'Tol', TOL);
check_norm('[alpha1; beta1]', '[alpha2; beta2]', 'Log', -1);

figure
subplot(121); res1.draw(PLOT_PROPS{:}); title('old')
subplot(122); res2.draw(PLOT_PROPS{:}); title('new')

%%
m = Mesh('../data/round_cuber.off');
V = m.V;
F = m.F;
% or: 
%mesh = Mesh(V, F);

H = m.H;
beta_g = m.generator_defects;

nv = m.nV; nf = m.nF;

rng(SEED)
[alpha1, beta1] = IOQ_highgenus_gpu(V, F);
res1 = TCODS(m, ...
             'k', [alpha1; beta1], ...
             'f0', FACE0, ...
             'theta0', THETA0, ...
             'degree', DEGREE, ...
             'CreateFField', true, ...
             'Duplicate', true, ...
             'gConstraintVec', GVEC);

rng(SEED)
[alpha2, beta2, x2, stats2, out2] = IOQ(V, F);
res2 = TCODS(m, ...
             'k', [alpha2; beta2], ...
             'f0', FACE0, ...
             'theta0', THETA0, ...
             'degree', DEGREE, ...
             'CreateFField', true, ...
             'Duplicate', true, ...
             'gConstraintVec', GVEC);

check_norm('res1.miq_energy', 'norm(x2)^2', 'Log', -1, 'Tol', TOL);
check_norm('[alpha1; beta1]', '[alpha2; beta2]', 'Log', -1);

figure
subplot(121); res1.draw(PLOT_PROPS{:}); title('old')
subplot(122); res2.draw(PLOT_PROPS{:}); title('new')
%set(gcf, 'WindowStyle', 'docked')