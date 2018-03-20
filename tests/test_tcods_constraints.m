FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = false;
fp = '../data/sphere_s0.off';
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;
CONSTRAINTS = [1, 0; 50, pi];

%% run IOQ
[alpha, beta, x, stats, out] = IOQ(...
    V, F, ...
    'GPU', USE_GPU, ...
    'Iterations', 2000, ...
    'Constraints', CONSTRAINTS);

disp(norm(x)^2)

%% run TCODS
k = [alpha; beta];

res_tc = TCODS(m, ...
    'k', k, ...
    'f0', FACE0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'Constraints', CONSTRAINTS);
%[E_edges, Emiq] = per_edge_energy(res_tc);

theta = res_tc.ffield_angles;
for i = 1:size(CONSTRAINTS, 2)
    fid = CONSTRAINTS(i, 1);
    t = CONSTRAINTS(i, 2);
    assert(norm(theta(fid) - t) < 1e-10)
end

disp(res_tc.miq_energy)

%%
figure
res_tc.draw('FaceAlpha', 0.8, 'FaceColor', 'w', 'EdgeColor', 'k')
hold on
res_tc.drawLabels()
hold off