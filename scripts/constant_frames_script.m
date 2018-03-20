%%
experiment_folder = 'constant_frames';
mkdir(experiment_folder);

p = find_data_folder();
FNAME = 'squares5_row.off';
file_path = fullfile(p, FNAME);
N = 4;
%constrained_faces = [1, 19];
%constraint_vectors = normalize_rows([0, 1, 0; 1, 1, 0]);
constrained_faces = [1];
constraint_vectors = normalize_rows([0, 1, 0]);

m = Mesh();
m.loadTM(file_path);
scale_factor = m.avg_length / 4;
                              
%% Direct rounding + PU
print_header('Direct rounding + PU')
solver = NRosy(m, N, constrained_faces, constraint_vectors);
solver.local_frames = [repmat([0, 1, 0], m.nF, 1); repmat([1, 0, 0], m.nF, 1)];
solver.k = zeros(m.nE, 1);
solver.create_NRosy_system();
solver.solve_with_direct_rounding();
solver.edgelist_phase_unwrapping2();

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

%% Direct rounding + PU
% print_header('Direct rounding')
% solver = NRosy(m, N, constrained_faces, constraint_vectors);
% solver.solve_with_direct_rounding();
% solver.edgelist_phase_unwrapping();
% 
% E = num2str(solver.E_quadratic);
% disp(['E : ', E]);
% disp('Plotting...')
% figure
% m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0)
% hold on
% solver.draw_constraints();
% solver.draw_field(scale_factor);
% title(['Direct Rounding + PU, E = ', num2str(solver.E_quadratic)])
% %m.drawLabels();
% hold off
% view([0,0,1])
% axis off
% path = fullfile(experiment_folder, ['5squares_row_direct_rounding_PU_E_', num2str(E)]);
% saveas(gcf, [path, '.png']);
% saveas(gcf, [path, '.fig']);