FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = false;

fp = '../data/round_cuber.off';
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;

CONSTRAINTS = [1, 0; 50, 0.1; 100, 0.2];
cfaces = CONSTRAINTS(:, 1);
cangs  = CONSTRAINTS(:, 2);
cvecs  = langles_to_vecs3d(m, cfaces, cangs);

[local_frames, frame_diffs] = create_local_frames(m);
[R, gamma_g] = create_constraints_mat(m, CONSTRAINTS(:, 1), CONSTRAINTS(:, 2), frame_diffs);
R = R';

PLOT_PROPS = {'FaceColor', 'w', 'PlotField', true, 'Constraints', CONSTRAINTS};

%% run IOQ
[alpha, beta, x, stats, out] = IOQ(...
    V, F, ...
    'GPU', USE_GPU, ...
    'Iterations', 2000, ...
    'Constraints', CONSTRAINTS);
gamma = out.gamma;
disp(norm(x)^2)
res_ioq = out.m;

%% run TCODS
k = [alpha; beta];

res_tc = TCODS(m, ...
    'k', k, ...
    'f0', FACE0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'Constraints', CONSTRAINTS, ...
    'Gamma', gamma);
%[E_edges, Emiq] = per_edge_energy(res_tc);

% TODO check with gamma
theta = res_tc.ffield_angles;
tangs = theta(cfaces);
for i = 2:size(CONSTRAINTS, 1)
    fid = CONSTRAINTS(i, 1);
    t = CONSTRAINTS(i, 2);
    %assert(norm(theta(fid) - (t + pi/2 * gamma(i-1))) < 1e-10)
    %fprintf('theta(%d) = %g\nt + pi/2*gamma(%d)=%g\n', fid, theta(fid), i-1, t+pi/2*gamma(i-1));
    check_norm('theta(fid)', 't+pi/2*gamma(i-1)');
end

disp(res_tc.miq_energy)

%%
[res_miq, elapsed_miq] = nrosy_mex(fp, cfaces, cvecs, DEGREE);

%%
title_ioq = sprintf('ioq, E = %g', res_ioq.miq_energy);
title_tc = sprintf('tc,  E = %g', res_tc.miq_energy);
title_miq = sprintf('miq, E = %g', res_miq.miq_energy);

figure
subplot(131); out.m.draw(PLOT_PROPS{:}); title(title_ioq)
subplot(132); res_tc.draw(PLOT_PROPS{:}); title(title_tc)
subplot(133); res_miq.draw(PLOT_PROPS{:}); title(title_miq)

%%
PLOT_PROPS = {'FaceColor', 'w', 'PlotField', true, 'Constraints', CONSTRAINTS, 'EdgeColor', 'k'};
figure
paths = {};
for i = 1:size(R, 2)
    paths{end+1} = find(R(:, i));
end
res_ioq.draw(PLOT_PROPS{:});
hold on; 
res_ioq.labelFaces(cfaces);
res_ioq.plotEdgePaths(paths, 'color', 'r'); 
hold off;
