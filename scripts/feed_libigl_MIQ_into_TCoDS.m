%% Settings
FNAME = 'sphere_s1.off';
OUT_FOLDER = 'sphere';
VERBOSE = true;
%I = [1, 1/4;    % Singularities
%     10, 1/4;
%     20, 1/4;
%     30, 1/4;
%     100, 1/4;
%     110, 1/4;
%     120, 1/4
%     130, 1/4];    
N = 4;          % Direction field degree (i.e., 4 for cross field)
f0 = [1];       % Starting face
v0 = [0, 0, 1]; % Starting direction
%v0 = [-0.8090,    0.3090,    0.5000];

% Find mesh file
mkdir(OUT_FOLDER);
p = find_data_folder();
%FNAME = 'sphere_s0.off';
%FNAME = 'torus_fat_r2.off';
fp = fullfile(p, FNAME);

% Load mesh
m = Mesh();
m.loadTM(fp);
scale = 4 * m.avg_length;

[d0, d1] = get_exterior_derivatives(m);
v0 = m.V(m.F(f0, 2), :) - m.V(m.F(f0, 1), :);
%v0 = [1, 1, 1];





%% libigl MIQ
print_header('libigl MIQ')
[R, S, periods, period_is_fixed, thetas, local_frames, k, EV, EF, FE] = nrosy_wrapper(fp, f0, v0, N);
solver_libigl = NRosy(m, N, f0, v0);
%solver_libigl.thetas = representative_to_local(R, solver_libigl.local_frames);
solver_libigl.k = k;
solver_libigl.thetas = thetas;
solver_libigl.periods = periods;
solver_libigl.local_frames = local_frames;
assert(norm(EF - m.EFAdj) == 0);
assert(norm(EV - m.EVAdj) == 0);
%assert(norm(FE - m.FEAdj) == 0);
figure(1)
subplot(224)
m.draw()
hold on
solver_libigl.draw_constraints(f0, v0);
solver_libigl.draw_field(scale, {'r', 'b', 'b', 'b'});
%draw_direction_field(m, R, find(abs(S)>1e-10), scale, {'r', 'b', 'b', 'b'});
title({'MIQ Greedy Rounding (libIGL)', ['E = ', num2str(solver_libigl.E_quadratic)]});
hold off


%%
print_header('MIQ direct rounding (my solver) with libigl periods')
solver = NRosy(m, N);
solver.period_is_fixed = ones(size(solver.period_is_fixed, 1), 1);
solver_libigl.k = k;
solver_libigl.thetas = thetas;
solver_libigl.periods = periods;
solver_libigl.local_frames = local_frames;
solver.create_NRosy_system(f0, v0, false);
solver.solve_with_direct_rounding();

subplot(221)
m.draw()
hold on
solver.draw_constraints(f0, v0);
solver.draw_field(scale, {'r', 'b', 'b', 'b'});
solver.E_quadratic();
title(['Direct Rounding w/ libigl p, E = ', num2str(solver.E_quadratic)])
hold off


%%
print_header('MIQ direct rounding (my solver)')
solver2 = NRosy(m, N, f0, v0);
solver2.create_NRosy_system(f0, v0, true);
solver2.solve_with_direct_rounding();
subplot(222)
m.draw()
hold on
solver2.draw_constraints(f0, v0);
solver2.draw_field(scale, {'r', 'b', 'b', 'b'});
title(['Direct rounding, E = ', num2str(solver2.E_quadratic)])
hold off








%% Compute trivial connection and build from it the direction field
inds = find(abs(S) > 1e-10);
I = [inds, -S(inds)];
connection = create_trivial_connection(m, I, N, VERBOSE);
DF = create_direction_field(m, connection, f0, v0, N);
subplot(223)
m.draw()
hold on
draw_direction_field(m, DF, I, scale, {'r', 'b', 'b', 'b'});
title('Trivial connection w/ libigl singularities')
hold off

%%
figure(2)
clf(2)
ax1 = subplot(121);
m.draw()
hold on
draw_direction_field(m, DF(0*m.nF+1:1*m.nF,:), I, scale, {'r', 'b', 'b', 'b'});
title('Trivial connection w/ libigl singularities')
hold off
ax2 = subplot(122);
m.draw()
hold on
draw_direction_field(m, R, find(abs(S)>1e-10), scale, {'r', 'b', 'b', 'b'});
title('libigl MIQ')
hold off
hlink = linkprop([ax1, ax2], {'CameraPosition', 'CameraUpVector'});

% 
% % MIQ Direct rounding
% %solver = NRosy(m, N, f0, v0);
% %solver.solve_with_direct_rounding();
% 
% % Plot results w/ trivial connection
% figure; 
% m.draw(ones(m.nV, 1), 'EdgeAlpha', 0.2, 'FaceAlpha', 0.0)
% %subplot(121)
% m.draw();
% hold on
% draw_direction_field(m, DF, I, scale, {'r', 'g', 'g', 'g'});
% hold off
% title('Trivial connection')
% 
% % Plot results MIQ direct rounding
% %subplot(122)
% %m.draw()
% %hold on
% %solver.draw_constraints(f0, v0);
% %solver.draw_field(scale, {'r', 'b', 'b', 'b'});
% %title(['Direct Rounding, E = ', num2str(solver.E_quadratic)])
% %hold off