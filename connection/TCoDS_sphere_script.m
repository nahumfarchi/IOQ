% Settings
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
I = [1, 2];
N = 4;          % Direction field degree (i.e., 4 for cross field)
f0 = 1;         % Starting face
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
scale = 4*m.avg_length;

% Compute trivial connection and build from it the direction field
connection = create_trivial_connection(m, I, N, VERBOSE);
DF = create_direction_field(m, connection, f0, v0, N);

% MIQ Direct rounding
%solver = NRosy(m, N, f0, v0);
%solver.solve_with_direct_rounding();

% Plot results w/ trivial connection
figure; 
m.draw(ones(m.nV, 1), 'EdgeAlpha', 0.2, 'FaceAlpha', 0.0)
%subplot(121)
m.draw();
hold on
draw_direction_field(m, DF, I, scale, {'r', 'b', 'b', 'b'});
hold off
title('Trivial connection')

% Plot results MIQ direct rounding
%subplot(122)
%m.draw()
%hold on
%solver.draw_constraints(f0, v0);
%solver.draw_field(scale, {'r', 'b', 'b', 'b'});
%title(['Direct Rounding, E = ', num2str(solver.E_quadratic)])
%hold off