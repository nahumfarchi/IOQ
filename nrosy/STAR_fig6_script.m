experiment_folder = 'STAR_fig6';
mkdir(experiment_folder);

p = find_data_folder();
FNAME = 'squares5_row.off';
file_path = fullfile(p, FNAME);
N = 4;
constrained_faces = [1, 19];
constraint_vectors = normalize_rows([0, 1, 0; 1, 1, 0]);

m = Mesh();
m.loadTM(file_path);
scale_factor = m.avg_length / 4;
                              
%% Direct rounding
print_header('Direct rounding')
solver = NRosy(m, N, constrained_faces, constraint_vectors);
solver.solve_with_direct_rounding();

E = num2str(solver.E_quadratic);
disp(['E : ', E]);
disp('Plotting...')
figure
m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0)
hold on
solver.draw_constraints();
solver.draw_field(scale_factor);
title(['Direct Rounding, E = ', num2str(solver.E_quadratic)])
%m.drawLabels();
hold off
view([0,0,1])
axis off
path = fullfile(experiment_folder, ['5squares_row_direct_rounding_E_', num2str(E)]);
saveas(gcf, [path, '.png']);
saveas(gcf, [path, '.fig']);

%% Direct rounding + PU
print_header('Direct rounding')
solver = NRosy(m, N, constrained_faces, constraint_vectors);
solver.solve_with_direct_rounding();
solver.edgelist_phase_unwrapping();

E = num2str(solver.E_quadratic);
disp(['E : ', E]);
disp('Plotting...')
figure
m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0)
hold on
solver.draw_constraints();
solver.draw_field(scale_factor);
title(['Direct Rounding + PU, E = ', num2str(solver.E_quadratic)])
%m.drawLabels();
hold off
view([0,0,1])
axis off
path = fullfile(experiment_folder, ['5squares_row_direct_rounding_PU_E_', num2str(E)]);
saveas(gcf, [path, '.png']);
saveas(gcf, [path, '.fig']);