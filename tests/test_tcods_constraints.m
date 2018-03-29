FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = false;

%fp = '../data/round_cuber.off';
%fp = '../data/rounded_cube_keenan_r_3k.off';
fp = '../data/fandisk.off';
m = Mesh(fp);
verts = m.V; faces = m.F; nv = m.nV; ne = m.nE; nf = m.nF;

%CONSTRAINTS = [1, 0; 50, 0.1; 100, 0.2];
%cfaces = CONSTRAINTS(:, 1);
%cangs  = CONSTRAINTS(:, 2);
%cvecs  = langles_to_vecs3d(m, cfaces, cangs);

% Constraints
mm = MESH(fp(1:end-4));
[V, D] = SHAPEOP.shape_operator_R(mm);
cfaces = find(D > 0.015);
% In 3d coordinates
%cvecs = [cfaces, V(cfaces, :)];
cvecs = V(cfaces, :);
% In angles with respect to the local frames
%cangs = [cfaces, vecs3d_to_langles(m, cfaces, V(cfaces, :))];
cangs = vecs3d_to_langles(m, cfaces, V(cfaces, :));

CONSTRAINTS = [cfaces, cangs];

[local_frames, frame_diffs] = create_local_frames(m);
[R, gamma_g] = create_constraints_mat(m, CONSTRAINTS(:, 1), CONSTRAINTS(:, 2), frame_diffs);
R = R';

%check_norm('langles_to_vecs3d(m, CONSTRAINTS(:, 1), CONSTRAINTS(:, 2:end))', ...
%    'cvecs(:, 2:end)');

PLOT_PROPS = {'FaceColor', 'w', 'PlotField', true, 'Constraints', CONSTRAINTS};

%% run IOQ
rng(SEED)
[alpha, beta, x, stats, out] = IOQ(...
    verts, faces, ...
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
% theta = res_tc.ffield_angles;
% tangs = theta(cfaces);
% for i = 2:size(CONSTRAINTS, 1)
%     fid = CONSTRAINTS(i, 1);
%     t = CONSTRAINTS(i, 2);
%     %assert(norm(theta(fid) - (t + pi/2 * gamma(i-1))) < 1e-10)
%     %fprintf('theta(%d) = %g\nt + pi/2*gamma(%d)=%g\n', fid, theta(fid), i-1, t+pi/2*gamma(i-1));
%     check_norm('theta(fid)', 't+pi/2*gamma(i-1)');
% end

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

%% Run ioq with init alpha from miq
inds = res_miq.vert_sing(:, 1);
alpha_miq = zeros(nv, 1); alpha_miq(inds) = res_miq.vert_sing(:, 2)*DEGREE;
rng(SEED)
[alpha2, beta2, x2, stats2, out2] = IOQ(...
    verts, faces, ...
    'GPU', USE_GPU, ...
    'Iterations', 2000, ...
    'Constraints', CONSTRAINTS, ...
    'init_alpha', alpha_miq);
gamma2 = out2.gamma;
res_ioq_init_alpha_miq = out2.m;
disp(norm(x2)^2)
res_ioq2 = out2.m;

%% Compare energy plots of ioq with initial:
%   1. random alpha
%   2. miq alpha
figure
subplot(131); 
plot(1:length(stats.Ehist), stats.Ehist, 1:length(stats2.Ehist), stats2.Ehist)
title('E'); legend('random initial alpha', 'miq alpha'); %ylim([0, 20])
subplot(132);
plot(1:length(stats.Ea_hist), stats.Ea_hist, 1:length(stats2.Ea_hist), stats2.Ea_hist)
title('Ea'); legend('random initial alpha', 'miq alpha'); %ylim([0, 20])
subplot(133);
plot(1:length(stats.Ec_hist), stats.Ec_hist, 1:length(stats2.Ec_hist), stats2.Ec_hist)
title('Ec'); legend('random initial alpha', 'miq alpha'); %ylim([0, 20])

%% Compare ioq w/ init miq alpha vs. miq
title_ioq_init_alpha_miq = sprintf('ioq w/ init miq alpha, E = %g', res_ioq_init_alpha_miq.miq_energy);
figure
subplot(131); res_ioq.draw(PLOT_PROPS{:}); title(title_ioq)
subplot(132); res_ioq_init_alpha_miq.draw(PLOT_PROPS{:}); title(title_ioq_init_alpha_miq)
subplot(133); res_miq.draw(PLOT_PROPS{:}); title(title_miq)

%%