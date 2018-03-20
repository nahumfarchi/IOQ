% experiment_folder = 'STAR_fig6';
% mkdir(experiment_folder);
% 
% p = find_data_folder();
% N = 4;
% 
% %constrained_faces = [1, 19];
% %constrained_faces = [1, 15];
% 
% %FNAME = 'squares1_row.off';
% %constrained_faces = [1, 3];
% %constraint_vectors = normalize_rows([0, 1, 0; 0, 1, 0]);
% %constraint_vectors = normalize_rows([0, 1, 0; 1, 1, 0]);
% 
% %FNAME = 'squares2_row.off';
% %constrained_faces = [1, 7];
% %constraint_vectors = normalize_rows([0, 1, 0; 1, 1, 0]);
% 
% FNAME = 'squares5_row.off';
% constrained_faces = [1, 19];
% constraint_vectors = normalize_rows([0, 1, 0; 1, 1, 0]);
% 
% file_path = fullfile(p, FNAME);
% 
% m = Mesh();
% m.loadTM(file_path);
% scale_factor = m.avg_length / 4;
% 
% d0 = create_d0(m);
% %H = create_holonomy(m); TODO
% A = [d0'];
% 
% % Plot figure with labels
% %subplot(222)
% m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0)
% m.drawLabels()
% view([0,0,1]); axis off;

%%
face0 = 1;
theta0 = 0;
degree = 1;
LOG = -1;

%% Genus 0
% m = Mesh('../data/bunny.off');
% 
% singularities = [10, -0.5; ...
%                  20, -0.5; ...
%                  30, 1.5; ...
%                  40, 0.5; ...
%                  50, 1];
% 
% 
% bunny_tcods = TCODS(m, singularities, [], face0, theta0, degree);
% bunny_tcods_mex = TCODS('../data/bunny.off', singularities, [], face0, theta0, degree, 'Mex', true);
% 
% figure(); 
% subplot(121); bunny_tcods.draw(); title(num2str(bunny_tcods.miq_energy));
% subplot(122); bunny_tcods_mex.draw(); title(num2str(bunny_tcods_mex.miq_energy));

%% Genus 1
fp = '../data/torus_s0.off';
m = Mesh(fp);
disp('generator_angle_defects: ')
generator_angle_defects(m)
disp('flatten_generators: ')
flatten_generators(m)
vert_sing = [1, 0;
             100, 0];
gen_sing = [1, 0; ...
            2, 0];
torus_tcods = TCODS(m, vert_sing, gen_sing, face0, theta0, degree);
x1 = torus_tcods.connection;
torus_tcods_mex = TCODS(fp, vert_sing, gen_sing, ...
    face0, theta0, degree, 'Mex', true);
x2 = torus_tcods.connection;
torus_tcods_mex.ffield_vectors = torus_tcods_mex.ffield_vectors ./ row_norm(torus_tcods_mex.ffield_vectors);
torus_tcods_mex.ffield_vectors = [FF(:, 3), FF(:, 2), FF(:, 1)];

figure();
subplot(121);
torus_tcods.draw();
title(['tcods ', num2str(torus_tcods.miq_energy)]);
subplot(122); 
torus_tcods_mex.draw(); 
title(['tcods keenan ', num2str(torus_tcods_mex.miq_energy)]);

% plot generator 1-forms
cycles = m.generator_cycles();
for i = 1:numel(cycles)
    cy = cycles{i};
    inds = find(m.H(:, i));
    edge_labels = m.H(inds, i);
    figure
    m.draw('FaceAlpha', 0.8, 'PlotField', true)
    hold on
    m.colorFaces(cy, 'g')
    m.labelEdges(inds, edge_labels, '%d', [0, 0, 0], 'color', 'r')
    m.labelEdges(inds, torus_tcods.connection(inds), '%.3g', [0, 0, 0.01], 'color', 'm')
    m.drawLabels('Edge', false, 'FontSize', 8)
    hold off
end

% plot flattened generators
[gen_Ad, flat_meshes] = flatten_generators(m);
[local_frames, frame_diffs] = create_local_frames(m);
for i = 1:numel(flat_meshes)
    mf = flat_meshes{i};
    [lf, fd] = create_local_frames(mf);
    inds = cycles{i};
    angles = torus_tcods.ffield_angles(inds);
    [ffield] = angles_to_ffield(angles, lf, degree);
    
    % Find which edges are participating in the generator cycle
    gen_edges = unique(reshape(m.FEAdj(inds, :), [], 1));
    %edge_labels = torus_tcods.connection(gen_edges);
    edge_labels = torus_tcods.H(gen_edges, i);
    
    figure()
    
    % df
    ax(1) = subplot(121);
    hold on
    mf.drawFaceField(ffield, 'AutoScaleFactor', 0.5*mf.avg_length)
    %mf.drawLabels('FontSize', 7)
    %mf.labelEdges(1:mf.nE-1, edge_labels, '%.3g', [0, 0, 0], 'color', 'b')
    %mf.labelEdges(1:mf.nE-1, frame_diffs(gen_edges), '%.3g', [0, 0.01, 0], 'color', 'r')
    mf.draw('FaceAlpha', 0.2)
    title('field')
    hold off
    
    % frames
    ax(2) = subplot(122);
    hold on
    mf.drawFaceField(lf(1:mf.nF, :), 'color', 'r', 'AutoScaleFactor', 0.5*mf.avg_length)
    mf.drawFaceField(lf(mf.nF+1:end, :), 'color', 'b', 'AutoScaleFactor', 0.5*mf.avg_length)
    %mf.drawLabels('FontSize', 7)
    mf.draw('FaceAlpha', 0.2)
    title('frames')
    hold off
    
    linkaxes(ax)
end

[d0, d1] = get_exterior_derivatives(m);
x = torus_tcods.connection;
check_norm('d1*x', '0', 'Log', LOG);

A = [d0'; m.H'];
check_norm('rank(full([A; d1]))', 'min(size([A; d1]))', 'Log', LOG);


%% Read keenan's system
A = load('torus_s0_A.txt');
A = sparse(A(:, 1)+1, A(:, 2)+1, A(:, 3), m.nE, m.nV+2);
A = A';
b = load('torus_s0_b.txt');
n = size(A, 2);
C = speye(n, n);
d = zeros(n, 1);
Aeq = A;
beq = b;
x_star = lsqlin(C, d, [], [], Aeq, beq, -inf(n, 1), inf(n, 1));
check_norm('A*x_star', 'b');
norm(x_star)^2

% QR decomposition
% A' E = Q R
% A x = -b
% E R' Q x = -b
% x = - (E R' Q)^{-1} b
[Q, R, E] = qr(A');
x_star = Q * (R' \ (E' * -b));
check_norm('A*x_star', '-b');
norm(x_star)^2

Kg = get_gaussian_curvature(m);
check_norm('b(1:m.nV)', '-Kg', 'Log', LOG);

% cycles
figure()
hold on

m.draw('FaceAlpha', 0.9, 'PlotField', false)

%m.plotEdgePaths({find(A(end-1,:)), find(A(end,:))})
path1 = find(A(end-1,:));
path2 = find(A(end,:));
%m.labelEdges(path1, A(end-1,path1))
m.labelEdges(path2, A(end,path2))
m.drawLabels('Edge', false)
hold off

%%
m = Mesh('../data/eight.off');

% cycles = m.generator_cycles;
% for i = 1:numel(cycles)
%     figure()
%     m.draw('FaceAlpha', 0.8, 'PlotField', false)
%     hold on
%     m.colorFaces(cycles{i}, 'g')
%     hold off
% end

disp('generator_angle_defects: ')
generator_angle_defects(m)
disp('flatten_generators: ')
flatten_generators(m, false)
vert_sing = [10, -1; ...
             100, -1];
gen_sing = [1, 0; ...
            2, 0; ...
            3, 0; ...
            4, 0];
eight_tcods = TCODS(m, vert_sing, gen_sing, face0, theta0, degree);
eight_tcods_mex = TCODS('../data/eight.off', vert_sing, gen_sing, face0, theta0, degree, 'Mex', true);
eight_tcods_mex.ffield_vectors = eight_tcods_mex.ffield_vectors ./ row_norm(eight_tcods_mex.ffield_vectors);
figure(); 
subplot(121); 
eight_tcods.draw(); 
title(['tcods ', num2str(eight_tcods.miq_energy)]);
subplot(122); 
eight_tcods_mex.draw(); 
title(['tcods keenan ', num2str(eight_tcods_mex.miq_energy)]);



% plot flattened generators
cycles = m.generator_cycles();
[gen_Ad, flat_meshes] = flatten_generators(m);
[local_frames, frame_diffs] = create_local_frames(m);
for i = 1:numel(flat_meshes)
    mf = flat_meshes{i};
    [lf, fd] = create_local_frames(mf);
    inds = cycles{i};
    angles = eight_tcods.ffield_angles(inds);
    [ffield] = angles_to_ffield(angles, lf, degree);
    
    % Find which edges are participating in the generator cycle
    gen_edges = unique(reshape(m.FEAdj(inds, :), [], 1));
    %edge_labels = torus_tcods.connection(gen_edges);
    edge_labels = eight_tcods.H(gen_edges, i);
    
    figure()
    
    % df
    ax(1) = subplot(121);
    hold on
    mf.drawFaceField(ffield, 'AutoScaleFactor', 2*mf.avg_length)
    %mf.drawLabels('FontSize', 7)
    %mf.labelEdges(1:mf.nE-1, edge_labels, '%.3g', [0, 0, 0], 'color', 'b')
    %mf.labelEdges(1:mf.nE-1, frame_diffs(gen_edges), '%.3g', [0, 0.01, 0], 'color', 'r')
    mf.draw('FaceAlpha', 0.2)
    title('field')
    hold off
    
    % frames
    ax(2) = subplot(122);
    hold on
    mf.drawFaceField(lf(1:mf.nF, :), 'color', 'r', 'AutoScaleFactor', 2*mf.avg_length)
    mf.drawFaceField(lf(mf.nF+1:end, :), 'color', 'b', 'AutoScaleFactor', 2*mf.avg_length)
    %mf.drawLabels('FontSize', 7)
    mf.draw('FaceAlpha', 0.2)
    title('frames')
    hold off
    
    linkaxes(ax)
end



[d0, d1] = get_exterior_derivatives(m);
x = eight_tcods.connection;
check_norm('d1*x', '0', 'Log', LOG);

A = [d0'; m.H'];
check_norm('rank(full([A; d1]))', 'min(size([A; d1]))', 'Log', LOG);