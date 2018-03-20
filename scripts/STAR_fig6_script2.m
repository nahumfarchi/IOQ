p = find_data_folder();
FNAME = 'squares5_row.off';
fp = fullfile(p, FNAME);

N = 4;
%f0 = [1, 3];
%f0 = [1, 7];
f0 = [1, 19];
%v0 = normalize_rows([0, 1, 0; 1, 1, 0]);
tc1 = pi+0.1;
tc2 = pi-0.1;
v0 = [cos(tc1), sin(tc1), 0; cos(tc2), sin(tc2), 0];

m = Mesh();
m.loadTM(fp);
scale_factor = m.avg_length / 6;

%% (1) libigl MIQ greedy rounding
title1 = 'libigl MIQ greedy rounding';
print_header(title1)
[R1, S1, periods1, period_is_fixed1, thetas1, local_frames1, r1, EV1, EF1, FE1] = nrosy_wrapper(fp, f0, v0, N);
E1 = E_MIQ(m, thetas1, r1, periods1, N);
DF1 = angles_to_DF(thetas1, local_frames1, N);
disp(['E1 = ', num2str(E1)])

[local_frames2, r2] = create_local_frames(m);

%% (2) MIQ direct rounding
title2 = 'MIQ Direct Rounding';
local_frames2 = [repmat([1, 0, 0], m.nF, 1); repmat([0, 1, 0], m.nF, 1)];
r2 = zeros(m.nE, 1);
print_header(title2)
[A2, b2, thetas2, periods2, theta_tags2, period_tags2] = create_MIQ_system(m, N, f0, v0, local_frames2, r2); 
[thetas2, periods2, E2] = solve_MIQ(A2, b2, thetas2, periods2, theta_tags2, period_tags2);
[thetas2_real, periods2_real, E2_real] = solve_MIQ(A2, b2, thetas2, periods2, theta_tags2, period_tags2, false);
DF2 = angles_to_DF(thetas2, local_frames2, N);
disp(['E2 = ', num2str(E2)])

%% (3) Greedy rounding + PU
title3 = 'Greedy rounding + PU';
print_header(title3)
grid_V = (1:m.nF)';
grid_E = m.EFAdj;
wrapped = thetas1;
%[unwrapped, n, P, Q] = unwrap_phases_grid(grid_V, grid_E, wrapped, N);
[unwrapped, n, edge_periods] = unwrap_phases(m, wrapped, r1, N);
DF3 = angles_to_DF(unwrapped, local_frames1, N);
E3 = E_MIQ(m, unwrapped, r1, periods1, N);
disp(['E3 = ', num2str(E3)])



%% Plots
figure(1); clf(1)
cf_colors = {'k', 'k', 'k', 'k'};
face_alpha = 0.0;
edge_alpha = 1.0;

subplot(223)
hold on
draw_constraints(m, f0, v0);
draw_direction_field(m, DF1, [], scale_factor, cf_colors, face_alpha, edge_alpha, periods1)
title({title1, ['E = ', num2str(E1)]})
view([0,0,1]); axis off;
hold off

subplot(222)
hold on
draw_constraints(m, f0, v0);
draw_direction_field(m, DF2, [], scale_factor, cf_colors, face_alpha, edge_alpha, periods2)
title({title2, ['E = ', num2str(E2)]})
view([0,0,1]); axis off;
hold off

subplot(221)
hold on
draw_constraints(m, f0, v0);
draw_direction_field(m, DF2, [], scale_factor, cf_colors, face_alpha, edge_alpha, periods2_real)
title({['MIQ Before Direct Rounding'], ['E = ', num2str(E2_real)]})
view([0,0,1]); axis off;
hold off

subplot(224)
hold on
m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0)
draw_local_frames(m, local_frames2)
m.draw(ones(m.nV, 1), 'FaceAlpha', 0.0, 'EdgeColor', period_tags2 > -1)
title('Local Frames'); hold off; view([0,0,1]); axis off
hold off

figure(2)
m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0)
m.drawLabels()
view([0,0,1]); axis off;