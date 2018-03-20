% Compare the results of
% 1) Direct rounding
% 2) Direct rounding + phase unwrapping
% 3) MIQ
% 4) MIQ + phase unwrapping

experiment_folder = 'direct_rounding_PU_direct_rounding_thetas_fixed';
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

%% Direct rounding
tic
print_header('Direct rounding')
solver = NRosy(m, N, constrained_faces, constraint_vectors);
solver.solve_with_direct_rounding();
toc

E = num2str(solver.E_quadratic);
disp(['E : ', E]);
disp('Plotting...')
figure
m.draw(ones(m.nV, 1), 'EdgeAlpha', 0.0)
hold on
solver.draw_constraints(constrained_faces, constraint_vectors);
solver.draw_field();
title(['Direct Rounding, E = ', num2str(solver.E_quadratic)])
%m.drawLabels();
hold off
path = fullfile(experiment_folder, ['bumpy_direct_rounding_E_', num2str(E)]);
saveas(gcf, [path, '.png']);
saveas(gcf, [path, '.fig']);
path = fullfile(experiment_folder, ['bumpy_direct_rounding_E_', num2str(E), '.rosy']);
solver.save(path);
%export_fig(path, '-png');
%export_fig(path, '-fig');

dbl_newline();

%% Direct rounding + PU
tic
print_header('Direct rounding + PU')
solver = NRosy(m, N, constrained_faces, constraint_vectors);
solver.solve_with_direct_rounding();
solver.edgelist_phase_unwrapping2();
toc

E = num2str(solver.E_quadratic);
disp(['E : ', E]);
disp('Plotting...')
figure
m.draw(ones(m.nV, 1), 'EdgeAlpha', 0.0)
hold on
solver.draw_constraints(constrained_faces, constraint_vectors);
solver.draw_field();
title(['Direct Rounding + Phase Unwrapping, E = ', num2str(solver.E_quadratic)])
hold off
path = fullfile(experiment_folder, ['bumpy_direct_rounding_PU_E_', num2str(E)]);
saveas(gcf, [path, '.png']);
saveas(gcf, [path, '.fig']);
path = fullfile(experiment_folder, ['bumpy_direct_rounding_PU_E_', num2str(E), '.rosy']);
solver.save(path);

dbl_newline();

%% MIQ (using IGL output)
print_header('IGL MIQ')
solver.thetas = importdata(fullfile(p, 'bumpy_igl_angles.txt'));
solver.periods = importdata(fullfile(p, 'bumpy_igl_periods.txt'));
E = num2str(solver.E_quadratic);
disp(['E : ', E]);
figure
m.draw(ones(m.nV, 1), 'EdgeAlpha', 0.0)
hold on
solver.draw_constraints(constrained_faces, constraint_vectors);
solver.draw_field();
title(['IGL, E = ', num2str(solver.E_quadratic)]);
hold off
path = fullfile(experiment_folder, ['bumpy_IGL_MIQ_E_', num2str(E)]);
saveas(gcf, [path, '.png']);
saveas(gcf, [path, '.fig']);
path = fullfile(experiment_folder, ['bumpy_IGL_MIQ_E_', num2str(E), '.rosy']);
solver.save(path);

dbl_newline();

%% MIQ + PU
print_header('IGL MIQ + PU')
tic
solver.edgelist_phase_unwrapping();
toc
E = num2str(solver.E_quadratic);
disp(['E : ', E]);

figure
m.draw(ones(m.nV, 1), 'EdgeAlpha', 0.0)
hold on
solver.draw_constraints(constrained_faces, constraint_vectors);
solver.draw_field();
title(['IGL MIQ + Phase Unwrapping', ', E = ', num2str(solver.E_quadratic)])
hold off
path = fullfile(experiment_folder, ['bumpy_IGL_MIQ_PU_E_', num2str(E)]);
saveas(gcf, [path, '.png']);
saveas(gcf, [path, '.fig']);
path = fullfile(experiment_folder, ['bumpy_IGL_MIQ_PU_E_', num2str(E), '.rosy']);
solver.save(path);

dbl_newline();