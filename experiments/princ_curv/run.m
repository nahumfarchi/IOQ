FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = false;
%fp = '../data/bunny.off';
%fp = '../data/sphere_s0.off';
fp = '../../../data/fandisk.off';
m = Mesh(fp);
mm = MESH(fp(1:end-4));
verts = m.V; faces = m.F; nv = m.nV; ne = m.nE; nf = m.nF;
%CONSTRAINTS = [1, 0; 50, 0.1; 100, 0.2; 150, 0.3; 110, 0.4; 120, 0.5];
%CONSTRAINTS = [1, 0];
%CONSTRAINTS = [1, 0; 50, 0.1];
[local_frames, frame_diffs] = create_local_frames(m);

% Constraints
[V, D] = SHAPEOP.shape_operator_R(mm);
inds = find(D > 0.2);
% In 3d coordinates
c3d = [inds, V(inds, :)];
% In angles with respect to the local frames
clangles = [inds, vecs3d_to_langles(m, inds, V(inds, :))];

check_norm('langles_to_vecs3d(m, clangles(:, 1), clangles(:, 2:end))', ...
    'c3d(:, 2:end)');

[R, gamma_g] = create_constraints_mat(m, clangles(:, 1), clangles(:, 2), frame_diffs);
R = R';

%% run IOQ
rng(SEED)
[alpha1, beta1, x1, stats, out] = IOQ(...
    verts, faces, ...
    'GPU', USE_GPU, ...
    'Iterations', 2000, ...
    'Constraints', clangles);
gamma1 = out.gamma;
res_ioq = out.m;

% if size(clangles, 1) > 1
%     for i = 2:size(clangles, 2)
%         fid = clangles(i, 1);
%         t = clangles(i, 2);
%         assert(norm(out.theta(fid) - (t + (pi/2)*out.gamma(i-1))) < 1e-8)
%     end
% end

disp(['Eioq = ', num2str(norm(x1)^2)])

%% run TCODS
k = [alpha1; beta1];

res_tc = TCODS(m, ...
    'k', k, ...
    'f0', FACE0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'Constraints', clangles, ...
    'Gamma', gamma1);
%[E_edges, Emiq] = per_edge_energy(res_tc);

% theta = res_tc.ffield_angles;
% for i = 1:size(clangles, 1)
%     fid = clangles(i, 1);
%     t = clangles(i, 2);
%     assert(norm(theta(fid) - t) < 1e-10)
% end

disp(['Etc = ', num2str(res_tc.miq_energy)])

%% Run MIQ
cfaces = c3d(:, 1);
cvecs  = c3d(:, 2:end);
[res_miq, elapsed_miq] = ...
    nrosy_mex(fp, cfaces, cvecs, DEGREE);

%% Plot energies
figure
hold on
plot(1:length(stats.Ea_hist), stats.Ea_hist)
plot(1:length(stats.Ec_hist), stats.Ec_hist)
plot(1:length(stats.Ehist),   stats.Ehist)
legend('Ea', 'Ec', 'E')

%% Compare hodge decomps
% ioq
[d0, d1] = get_exterior_derivatives(m);
d0d0f = factorize(d0'*d0);
d1d1f = factorize(d1*d1');
a1 = out.a; b1 = out.b; c1 = out.c;
a11 = d0d0f \ (d0' * x1);
c11 = d1d1f \  (d1 * x1);
check_norm('d0*a1', 'd0*a11', 'Log', -1);
check_norm('d1''*c1', 'd1''*c11', 'Log', -1);

Ea1 = norm(d0*a1)^2;
Ec1 = norm(d1'*c1)^2;
E1 = Ea1 + Ec1;
check_norm('E1', 'res_ioq.miq_energy', 'Log', -1);

% tc
x2 = res_tc.connection;
a2 = d0d0f \ (d0' * x2);
c2 = d1d1f \ (d1 * x2);

Ea2 = norm(d0*a2)^2;
Ec2 = norm(d1'*c2)^2;
E2 = Ea2 + Ec2;
check_norm('E2', 'res_tc.miq_energy', 'Log', -1);

% miq
theta3 = res_miq.ffield_angles;
p3 = res_miq.periods;
r3 = res_miq.frame_diffs;

% change of vars (theta3, p3, r3) --> (x3, alpha3, beta3)
%alpha3 = round((2/pi) * (alpha_g - d0'*r3) - d0'*p3);
%beta3 = round((2/pi) * (beta_g - H'*r3) - H'*p3);
x3 = -r3 - d1'*theta3 - (pi/2)*p3;
a3 = d0d0f \ (d0' * x3);
c3 = d1d1f \ (d1 * x3);

Ea3 = norm(d0*a3)^2;
Ec3 = norm(d1'*c3)^2;
E3 = Ea3 + Ec3;
check_norm('E3', 'res_miq.miq_energy');

% print info
fprintf('ioq\n----\nE = %g\nEa = %g\nEc = %g\n\n', E1, Ea1, Ec1);
fprintf('tc\n-----\nE = %g\nEa = %g\nEc = %g\n\n', E2, Ea2, Ec2);
fprintf('miq\n-----\nE = %g\nEa = %g\nEc = %g\n', E3, Ea3, Ec3);

% tests
check_norm('R''*x1', '(pi/2)*gamma1 - gamma_g', 'Log', -1);
check_norm('R''*x2', '(pi/2)*gamma1 - gamma_g', 'Log', -1);

%% Plot constraint paths

%[R, gamma_g] = create_constraints_mat(m, clangles(:, 1), clangles(:, 2), frame_diffs);
%R = R';
paths = {};
for i = 1:size(R, 2)
    paths{end+1} = find(R(:, i));
end

%%
figure
func = zeros(nf, 1);
func(c3d(:, 1)) = 1;
res_ioq.draw('Constraints', [], 'PlotSing', false, 'Func', func, 'EdgeColor', 'k');
hold on; res_ioq.plotEdgePaths(paths, 'color', 'r'); hold off;
%subplot(122)
%res_ioq.draw('FaceColor', 'w', 'Constraints', c3d, 'EdgeColor', 'k', 'PlotSing', false);

%%

%plot_props = {'FaceColor', 'flat', 'EdgeColor', 'k', 'PlotField', false, ...
%    'Func', mm.Ie2f*x1, ...
%    'Constraints', c3d};
plot_props = {'FaceColor', 'flat', 'PlotField', true, 'PlotSing', true, ...
    'nrosy_colors', {'w', 'w', 'w', 'w'}, 'Colorbar', false, ...
    'Constraints', [], 'MarkerSize', 20};

title_ioq = {'ioq', sprintf('E = %g, Ea = %g, Ec = %g, ns = %d', res_ioq.miq_energy, Ea1, Ec1, res_tc.n_singularities)};
title_tc = {'tc', sprintf('E = %g, Ea = %g, Ec = %g, ns = %d', res_tc.miq_energy, Ea2, Ec2, res_ioq.n_singularities)};
title_miq = {'miq', sprintf('E = %g, Ea = %g, Ec = %g, ns = %d', res_miq.miq_energy, Ea3, Ec3, res_miq.n_singularities)};

figure
subplot(131); res_ioq.draw(plot_props{:}); title(title_ioq); 
subplot(132); res_tc.draw(plot_props{:}); title(title_tc); 
subplot(133); res_miq.draw(plot_props{:}); title(title_miq);
%hold on
%res_tc.drawLabels()
%hold off

%%
plot_props = {'FaceColor', 'flat', 'PlotField', true, 'PlotSing', true, ...
    'nrosy_colors', {'w', 'w', 'w', 'w'}, 'Colorbar', false, ...
    'Constraints', [], 'MarkerSize', 20};
plot_props2 = {'Func', mm.Ie2f*x1.^2, plot_props{:}};
[E, E_edges] = E_MIQ(m, res_miq.ffield_angles, res_miq.frame_diffs, res_miq.periods, DEGREE);
plot_props3 = {'Func', mm.Ie2f*E_edges, plot_props{:}};

figure
ax(1) = subplot(121); res_ioq.draw(plot_props2{:}); title(title_ioq); 
hold on; res_ioq.plotEdgePaths(paths, 'color', 'r'); hold off;
ax(2) = subplot(122); res_miq.draw(plot_props3{:}); title(title_miq); 
linkaxes(ax)

return


%% Compute c using our method and the MIQ singularities
rng(SEED)
alpha4 = zeros(nv, 1); alpha4(res_miq.vert_sing(:, 1)) = DEGREE*res_miq.vert_sing(:, 2);
[alpha4, beta4, x4, stats4, out4] = IOQ(...
    verts, faces, ...
    'GPU', USE_GPU, ...
    'Iterations', 2000, ...
    'Constraints', clangles, ...
    'alpha', alpha4, ...
    'GatherStats', false);
gamma4 = out4.gamma;
res_ioq4 = out4.m;
              
% Use miq a and c (hodge decomp)
gamma5 = ( (R' * d1')*c3 + gamma_g + R'*d0*a3 ) * 2 / pi;
[alpha5, beta5, x5, stats5, out5] = IOQ(...
    verts, faces, ...
    'GPU', USE_GPU, ...
    'Iterations', 2000, ...
    'Constraints', clangles, ...
    'alpha', alpha4, ...
    'gamma', gamma5, ...
    'GatherStats', false);
gamma5 = out5.gamma;
res_ioq5 = out5.m;
              
%%
title4 = sprintf('ioq, E = %g, ns = %d', res_ioq.miq_energy, res_ioq.n_vert_sing);
title5 = sprintf('ioq w/ miq vert sing, E = %g, ns = %d', res_ioq4.miq_energy, res_ioq4.n_vert_sing);
title6 = {'ioq w/ miq vert and const sing', sprintf('E = %g, ns = %d', res_ioq5.miq_energy, res_ioq5.n_vert_sing)};

figure
subplot(131); res_ioq.draw(plot_props{:}); title(title4); view(3)
subplot(132); res_ioq4.draw(plot_props{:}); title(title5); view(3)
subplot(133); res_ioq5.draw(plot_props{:}); title(title6); view(3)

%%
title7 = sprintf('miq, E = %g', res_miq.miq_energy);
figure;
subplot(121); res_miq.draw(plot_props{:}); title(title7); colorbar off
subplot(122); res_ioq4.draw(plot_props{:}); title(title5); colorbar off
%subplot(122); res_ioq5.draw(plot_props{:}); title(title6); colorbar off

%%
PLOT_PROPS = {'FaceColor', 'w', 'PlotField', true, 'Constraints', c3d, 'EdgeColor', 'k'};
figure
res_ioq.draw(PLOT_PROPS{:});
%%
paths = {};
for i = 1:size(R, 2)
    paths{end+1} = find(R(:, i));
end
hold on; 
res_ioq.labelFaces(cfaces);
res_ioq.plotEdgePaths(paths, 'color', 'r'); 
hold off;





