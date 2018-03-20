% Compare the results of
% 1) Direct rounding
% 2) Direct rounding + phase unwrapping
% 3) MIQ
% 4) MIQ + phase unwrapping

experiment_folder = 'tmp3';
mkdir(experiment_folder);

p = find_data_folder();
screen_size = get(0, 'ScreenSize');
width = screen_size(3);
height = screen_size(4);

m = Mesh();
m.loadTM(fullfile(p, 'bumpy.off'));
N = 4;
% constrained_faces = [141, 1574];
% constraint_vectors = [ 1, 0, 0;
%                       0, 0, -1;];
constrained_faces = [1];
constraint_vectors = [1, 1, 1];

dbl_newline();

figure(FIG('MIQ_PU_script'));

%% Direct rounding
tic
print_header('Direct rounding')
solver = NRosy(m, N, constrained_faces, constraint_vectors);
solver.solve_with_direct_rounding();
toc

E = num2str(solver.E_quadratic);
disp(['E : ', E]);
disp('Plotting...')
subplot(221)
m.draw(ones(m.nV, 1), 'EdgeAlpha', 0.0)
hold on
solver.draw_constraints(constrained_faces, constraint_vectors);
solver.draw_field();
title({'MIQ Direct Rounding', ['E = ', num2str(solver.E_quadratic)]})
%m.drawLabels();
hold off
path = fullfile(experiment_folder, ['bumpy_direct_rounding_E_', num2str(E)]);
saveas(gcf, [path, '.png']);
saveas(gcf, [path, '.fig']);
%export_fig(path, '-png');
%export_fig(path, '-fig');

dbl_newline();

%% Direct rounding + PU
tic
print_header('Direct rounding + PU')
solver2 = NRosy(m, N, constrained_faces, constraint_vectors);
solver2.solve_with_direct_rounding();
solver2.edgelist_phase_unwrapping2();
toc

E = num2str(solver2.E_quadratic);
disp(['E : ', E]);
disp('Plotting...')
subplot(222)
m.draw(ones(m.nV, 1), 'EdgeAlpha', 0.0)
hold on
solver2.draw_constraints(constrained_faces, constraint_vectors);
solver2.draw_field();
title({'MIQ Direct Rounding + Phase Unwrapping', ['E = ', num2str(solver2.E_quadratic)]})
hold off
path = fullfile(experiment_folder, ['bumpy_direct_rounding_PU_E_', num2str(E)]);
saveas(gcf, [path, '.png']);
saveas(gcf, [path, '.fig']);

dbl_newline();

%% MIQ (using IGL output)
print_header('IGL MIQ')
solver3 = NRosy(m, N, constrained_faces, constraint_vectors);
solver3.thetas = importdata(fullfile(p, 'bumpy_igl_angles.txt'));
solver3.periods = importdata(fullfile(p, 'bumpy_igl_periods.txt'));
E = num2str(solver3.E_quadratic);
disp(['E : ', E]);
subplot(223)
m.draw(ones(m.nV, 1), 'EdgeAlpha', 0.0)
hold on
solver3.draw_constraints(constrained_faces, constraint_vectors);
solver3.draw_field();
title({'MIQ Greedy Rounding (libIGL)', ['E = ', num2str(solver3.E_quadratic)]});
hold off
path = fullfile(experiment_folder, ['bumpy_IGL_MIQ_E_', num2str(E)]);
saveas(gcf, [path, '.png']);
saveas(gcf, [path, '.fig']);

dbl_newline();

%% MIQ + PU
print_header('IGL MIQ + PU')
tic
solver4 = NRosy(m, N, constrained_faces, constraint_vectors);
solver4.edgelist_phase_unwrapping2();
toc
E = num2str(solver4.E_quadratic);
disp(['E : ', E]);

subplot(224)
m.draw(ones(m.nV, 1), 'EdgeAlpha', 0.0)
hold on
solver4.draw_constraints(constrained_faces, constraint_vectors);
solver4.draw_field();
title({'MIQ Greedy Rounding (libIGL) + Phase Unwrapping', ['E = ', num2str(solver4.E_quadratic)]})
hold off
path = fullfile(experiment_folder, ['bumpy_IGL_MIQ_PU_E_', num2str(E)]);
saveas(gcf, [path, '.png']);
saveas(gcf, [path, '.fig']);

dbl_newline();