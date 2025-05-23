%% =======================================================================
%  Init IOQ with singularities from GO. Is convergence faster?
%  =======================================================================

FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = false;

fp = '../../../data/kitten_rg.off';
m = Mesh(fp);
V = m.V; F = m.F; ne = m.nE; nv = m.nV; nf = m.nF; ng = m.genus;

% Run GO
[res_go, Ego, elapsed_go] = GO(fp, 4);
alpha_go = zeros(nv, 1);
inds = res_go.vert_sing(:, 1);
alpha_go(inds) = res_go.vert_sing(:, 2);

% Run IOQ with random init singularities
rng(SEED);
[alpha, beta, elapsed1, ~, out1] = IOQ_highgenus_gpu(V, F, 'UseGPU', USE_GPU, 'highg_method', 'option1a');
k = [alpha; beta];
res1 = TCODS(m, ...
    'k', k, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);
%[E_edges, Emiq] = per_edge_energy(res);

% Run IOQ with GO singularities
rng(SEED);
[alpha, beta, elapsed2, ~, out2] = IOQ_highgenus_gpu(V, F, ...
    'UseGPU', USE_GPU, ...
    'highg_method', 'option1a', ...
    'alpha_P', alpha_go);
k = [alpha; beta];
res2 = TCODS(m, ...
    'k', k, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);