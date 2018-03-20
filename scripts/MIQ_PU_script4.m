%% Settings
%FNAME = 'sphere_s0.off';
FNAME = 'bumpy.off';
OUT_FOLDER = '21.5.17';
VERBOSE = true;

mkdir(OUT_FOLDER);

% Find and load mesh file
DATA_FOLDER = find_data_folder();
fp = fullfile(DATA_FOLDER, FNAME);

m = Mesh();
m.loadTM(fp);

scale = m.avg_length / 2;

N = 4;
f0 = [1];
v0 = m.V(m.F(f0, 2), :) - m.V(m.F(f0, 1), :);
v0 = v0 / norm(v0);

[d0, d1] = get_exterior_derivatives(m);


%% (1) libigl MIQ greedy rounding
title1 = '(1) libigl MIQ greedy rounding';
print_header(title1)
[R1, S1, periods1, period_is_fixed1, thetas1, local_frames1, r1, EV1, EF1, FE1] = nrosy_wrapper(fp, f0, v0, N);
E1 = E_MIQ(m, thetas1, r1, periods1, N);
DF1 = angles_to_DF(thetas1, local_frames1, N);
disp(['E1 = ', num2str(E1)])

[local_frames2, r2] = create_local_frames(m);


%% (2) MIQ direct rounding
title2 = '(2) MIQ Direct Rounding';
print_header(title2)
[A2, b2, thetas2, periods2, theta_tags2, period_tags2] = create_MIQ_system(m, N, f0, v0, local_frames1, r1); 
[thetas2, periods2, E2] = solve_MIQ(A2, b2, thetas2, periods2, theta_tags2, period_tags2);
DF2 = angles_to_DF(thetas2, local_frames1, N);
disp(['E2 = ', num2str(E2)])


%% (3) MIQ direct rounding + phase unwrapping
title3 = '(3) MIQ Direct Rounding + Phase Unwrapping';
print_header(title3)
[A3, b3, thetas3, periods3, theta_tags3, period_tags3] = create_MIQ_system(m, N, f0, v0, local_frames1, r1); 
[thetas3, periods3, E3] = solve_MIQ(A3, b3, thetas3, periods3, theta_tags3, period_tags3);
[thetas3_unwrapped, n, periods3_PU] = unwrap_phases(m, thetas2, r1, N);
E3 = E_MIQ(m, thetas3, r1, periods3_PU, N);
%E3 = E_MIQ(m, thetas3_unwrapped, r1, periods3, N);
DF3 = angles_to_DF(thetas3, local_frames1, N);
disp(['E3 = ', num2str(E3)])


%% Plots
figure(1); clf(1)
cf_colors = {'k', 'k', 'k', 'k'};

draw_direction_field(m, DF3, [], scale, cf_colors)