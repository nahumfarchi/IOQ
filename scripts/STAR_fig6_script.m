experiment_folder = 'STAR_fig6';
mkdir(experiment_folder);

p = find_data_folder();
N = 4;

%constrained_faces = [1, 19];
%constrained_faces = [1, 15];

%FNAME = 'squares1_row.off';
%constrained_faces = [1, 3];
%constraint_vectors = normalize_rows([0, 1, 0; 0, 1, 0]);
%constraint_vectors = normalize_rows([0, 1, 0; 1, 1, 0]);

%FNAME = 'squares2_row.off';
%constrained_faces = [1, 7];
%constraint_vectors = normalize_rows([0, 1, 0; 1, 1, 0]);

FNAME = 'squares5_row.off';
constrained_faces = [1, 19];
constraint_vectors = normalize_rows([0, 1, 0; 1, 1, 0]);

file_path = fullfile(p, FNAME);

m = Mesh();
m.loadTM(file_path);
scale_factor = m.avg_length / 4;
                              
%% Direct rounding
print_header('Direct rounding')
solver = NRosy(m, N, constrained_faces, constraint_vectors);

figure
subplot(221)
m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0)
hold on
solver.draw_constraints(constrained_faces, constraint_vectors);
solver.draw_field(scale_factor);
title(['Before, E = ', num2str(solver.E_quadratic)])
%m.drawLabels();
hold off
view([0,0,1])
axis off


solver.solve_with_direct_rounding();

E = num2str(solver.E_quadratic);
disp(['E : ', E]);
disp('Plotting...')


subplot(222)
m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0)
m.drawLabels()
view([0,0,1]); axis off;

subplot(223)
hold on
m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0)
draw_local_frames(m, solver.local_frames)
title('Local Frames'); hold off; view([0,0,1]); axis off

subplot(224)
m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0)
hold on
solver.draw_constraints(constrained_faces, constraint_vectors);
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
solver.draw_constraints(constrained_faces, constraint_vectors);
solver.draw_field(scale_factor);
title(['Direct Rounding + PU, E = ', num2str(solver.E_quadratic)])
%m.drawLabels();
hold off
view([0,0,1])
axis off
path = fullfile(experiment_folder, ['5squares_row_direct_rounding_PU_E_', num2str(E)]);
saveas(gcf, [path, '.png']);
saveas(gcf, [path, '.fig']);