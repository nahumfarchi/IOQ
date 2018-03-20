% Compare the results of
% 1) Direct rounding
% 2) Direct rounding + phase unwrapping
% 3) MIQ
% 4) MIQ + phase unwrapping

experiment_folder = 'direct_rounding_PU_iter';
mkdir(experiment_folder);

p = find_data_folder();

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
solver.draw_constraints();
solver.draw_field();
title(['Direct Rounding Rand Perm, E = ', num2str(solver.E_quadratic)])
%m.drawLabels();
hold off
path = fullfile(experiment_folder, ['bumpy_direct_rounding_E_', num2str(E)]);
saveas(gcf, [path, '.png']);
saveas(gcf, [path, '.fig']);
%export_fig(path, '-png');
%export_fig(path, '-fig');

dbl_newline();

%% Direct rounding + PU with iterations
tic
print_header('Direct rounding + PU')
solver = NRosy(m, N, constrained_faces, constraint_vectors);
solver.solve_with_direct_rounding();
solver.edgelist_phase_unwrapping();
toc

E = num2str(solver.E_quadratic);
disp(['E : ', E]);
disp('Plotting...')
figure
m.draw(ones(m.nV, 1), 'EdgeAlpha', 0.0)
hold on
solver.draw_constraints();
solver.draw_field();
title(['Rand Perm, Direct Rounding + PU, E = ', num2str(solver.E_quadratic)])
%m.drawLabels();
hold off
path = fullfile(experiment_folder, ['rand_perm_bumpy_direct_rounding_E_', num2str(E)]);
saveas(gcf, [path, '.png']);
saveas(gcf, [path, '.fig']);

dbl_newline();