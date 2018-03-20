VERBOSE = true;
FNAME = 'sphere_s1.off';
p = find_data_folder();
fp = fullfile(p, FNAME);


m = Mesh();
m.loadTM(fp);
scale = 4*m.avg_length;

N = 1;
local_frames = create_local_frames(m);
ks = [1, 2];                 % singularity indices
cfids_list = [1, 19];        % ids of constrained faces
cthetas_list = [pi/2, pi/4]; % constraint angles
cthetas = nan(m.nF, 1);
cthetas(cfids_list) = cthetas_list;
cvectors = local_angles_to_gvectors(m, cthetas, local_frames, N);
f0 = cfids_list(1);
theta0 = cthetas_list(1);


% Compute tc
connection = create_trivial_connection(m, ks, N, VERBOSE, cfids_list, cthetas_list);
%connection = create_trivial_connection(m, ks, N, VERBOSE);
DF = create_direction_field(m, connection, f0, theta0, N);

%% Plot
cf_colors = {'b', 'b', 'b', 'b'};
face_alpha = 1.0;
edge_alpha = 0.2;

figure(1); clf(1)
%m.draw(ones(m.nV, 1), 'EdgeAlpha', 0.5, 'FaceAlpha', 0.0)
hold on
draw_constraints(m, cfids_list, cvectors)
draw_direction_field(m, DF, ks, scale_factor, cf_colors, face_alpha, edge_alpha)
hold off
title('TC w/ multiple constraints')

