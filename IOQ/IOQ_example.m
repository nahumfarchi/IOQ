% =========================================================================
% IQO example
% =========================================================================
FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = false;
fp = '../../data/rounded_cube_keenan_r_3k.off';
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;

% run IOQ
[alpha, beta] = IOQ_highgenus_gpu(...
    V, F, ...
    'UseGPU', USE_GPU, ...
    'Iterations', 2000);

% run TCODS
k = [alpha; beta];
res_tc = TCODS(m, ...
    'k', k, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);
[E_edges, Emiq] = per_edge_energy(res_tc);

% Check that TC and MIQ energy are the same.
% The .miq_energy field is actually |x|^2
assert(abs(Emiq - res_tc.miq_energy) < 1e-10)

% Plot
figure
res_tc.draw('FaceAlpha', 1, 'EdgeAlpha', 0)
    