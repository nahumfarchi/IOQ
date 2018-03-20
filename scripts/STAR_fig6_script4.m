p = find_data_folder();
FNAME = 'sphere_s0.off';
%FNAME = 'sphere_s0.off';
fp = fullfile(p, FNAME);

N = 4;
%f0 = [1, 3];
%f0 = [1, 7];
fids_list = [1, 19];
%v0 = normalize_rows([0, 1, 0; 1, 1, 0]);
%tc1 = pi+0.1;
%tc2 = pi-0.1;
%tc1 = pi / 2 - 1e-10;
%tc2 = pi / 4;
tc1 = 0.75*pi;
tc2 = 0.1*pi;
thetas_list = [tc1, tc2];
%v0 = [cos(tc1), sin(tc1), 0; cos(tc2), sin(tc2), 0];
tmp_thetas = nan(m.nF, 1);
tmp_thetas(fids_list) = thetas_list;
v0 = local_angles_to_gvectors(m, tmp_thetas, local_frames2, 1);

m = Mesh();
m.loadTM(fp);
scale_factor = m.avg_length / 6;
[d0, d1] = get_exterior_derivatives(m);

%% (1) libigl MIQ greedy rounding
tic
title1 = 'libigl MIQ greedy rounding';
print_header(title1)
[R1, S1, periods1, period_is_fixed1, thetas1, local_frames1, r1, EV1, EF1, FE1] = nrosy_wrapper(fp, fids_list, v0, N);
E1 = E_MIQ(m, thetas1, r1, periods1, N);
DF1 = angles_to_DF(thetas1, local_frames1, N);
disp(['E1 = ', num2str(E1)])
toc

%% (2) MIQ direct rounding
title2 = 'MIQ Direct Rounding';
%local_frames2 = [repmat([1, 0, 0], m.nF, 1); repmat([0, 1, 0], m.nF, 1)];
%r2 = zeros(m.nE, 1);
[local_frames2, r2] = create_local_frames(m);
print_header(title2)
[A2, b2, thetas2, periods2, theta_tags2, period_tags2] = create_MIQ_system(m, N, fids_list, v0, local_frames2, r2); 
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

%% (4) min ||\eta||_1+||d0\theta+2\pi\eta||_2^2
%      s.t. \eta \in \mathbb{Z}
%           -\pi \leq \theta_i \leq \theta
%           \theta_i = \hat{\theta_i} \forall v_i \in V_c


%% (5) min ||x||_2^2 + ||d0^T*x
%      s.t. Cx = bc
tic
[C, bc] = create_constraints_mat(m, fids_list, thetas_list, r2);
lambda = 1.0;
cvx_begin
    variable x5(m.nE)
    minimize( norm(x5, 2) + lambda*pi/2*norm(d0'*x5+get_gaussian_curvature(m), 1) )
    %minimize( norm(x5, 2) + lambda*pi/2*norm(d0'*x5, 1) )
    subject to
        C*x5 == bc
cvx_end

[DF5, thetas5] = create_direction_field(m, x5, fids_list(1), thetas_list(1), N, local_frames2, r2);
E5 = cvx_optval;
title5 = ['|x|_2^2 + |d_0^Tx+G_k|_1 = ', num2str(norm(x5)), ' + ', num2str(norm(d0'*x5,1)), ' = ', num2str(E5)];
toc


%% Plots
ZOOM_FAC = 1.4;

figure(1); clf(1)
cf_colors = {'k', 'k', 'k', 'k'};
face_alpha = 0.0;
edge_alpha = 1.0;

subplot(425)
hold on
draw_constraints(m, fids_list, v0);
draw_direction_field(m, DF1, [], scale_factor, cf_colors, face_alpha, edge_alpha, periods1)
title({title1, ['|AX-b|_2^2 = ', num2str(E1, '%.4f'), ...
                ', |p|_2^2 = ', num2str(norm(periods1)^2, '%.4f')]})
view([0,0,1]); axis off;
hold off
zoom(ZOOM_FAC)
subplot(427)
m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0);
m.labelFaces(thetas1, '%.2f');
view([0,0,1]); axis off;
zoom(ZOOM_FAC)

subplot(422)
hold on
draw_constraints(m, fids_list, v0);
draw_direction_field(m, DF2, [], scale_factor, cf_colors, face_alpha, edge_alpha, periods2)
%title({title2, ['E = ', num2str(E2)]})
title({title2, ['|AX-b|_2^2 = ', num2str(E2, '%.4f'), ...
                ', |p|_2^2 = ', num2str(norm(periods2)^2, '%.4f')]})
view([0,0,1]); axis off;
hold off
zoom(ZOOM_FAC)
subplot(424)
m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0);
m.labelFaces(thetas2, '%.2f');
view([0,0,1]); axis off;
zoom(ZOOM_FAC)

subplot(421)
hold on
draw_constraints(m, fids_list, v0);
draw_direction_field(m, DF2, [], scale_factor, cf_colors, face_alpha, edge_alpha, periods2_real)
%title({['MIQ p\in\mathbb{R}'], ['E = ', num2str(E2_real)]})
title({'MIQ real p', ['|AX-b|_2^2 = ', num2str(E2_real, '%.4f'), ...
                ', |p|_2^2 = ', num2str(norm(periods2_real)^2, '%.4f')]})
view([0,0,1]); axis off;
hold off
zoom(ZOOM_FAC)
subplot(423)
m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0);
m.labelFaces(thetas2_real, '%.2f');
view([0,0,1]); axis off;
zoom(ZOOM_FAC)

subplot(426)
hold on
draw_constraints(m, fids_list, v0);
draw_direction_field(m, DF5, [], scale_factor, cf_colors, face_alpha, edge_alpha)
title(title5)
view([0,0,1]); axis off;
hold off
zoom(ZOOM_FAC)
subplot(428)
m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0);
m.labelFaces(thetas5, '%.2f');
view([0,0,1]); axis off;
zoom(ZOOM_FAC);

% subplot(248)
% hold on
% m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0)
% draw_local_frames(m, local_frames2)
% m.draw(ones(m.nV, 1), 'FaceAlpha', 0.0, 'EdgeColor', period_tags2 > -1)
% title('Local Frames'); hold off; view([0,0,1]); axis off
% hold off
% 
% figure(2)
% m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0)
% m.drawLabels()
% view([0,0,1]); axis off;

%% Plots
figure(3); clf(3)
cf_colors = {'k', 'k', 'k', 'k'};
face_alpha = 0.0;
edge_alpha = 1.0;

subplot(223)
hold on
draw_constraints(m, fids_list, v0);
draw_direction_field(m, DF1, [], scale_factor, cf_colors, face_alpha, edge_alpha, periods1)
title({title1, ['|AX-b|_2^2 = ', num2str(E1, '%.4f'), ...
                ', |p|_2^2 = ', num2str(norm(periods1)^2, '%.4f')]})
view([0,0,1]); axis off;
hold off
zoom(ZOOM_FAC)
% subplot(427)
% m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0);
% m.labelFaces(thetas1, '%.2f');
% view([0,0,1]); axis off;
% zoom(ZOOM_FAC)

subplot(222)
hold on
draw_constraints(m, fids_list, v0);
draw_direction_field(m, DF2, [], scale_factor, cf_colors, face_alpha, edge_alpha, periods2)
%title({title2, ['E = ', num2str(E2)]})
title({title2, ['|AX-b|_2^2 = ', num2str(E2, '%.4f'), ...
                ', |p|_2^2 = ', num2str(norm(periods2)^2, '%.4f')]})
view([0,0,1]); axis off;
hold off
zoom(ZOOM_FAC)
% subplot(424)
% m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0);
% m.labelFaces(thetas2, '%.2f');
% view([0,0,1]); axis off;
% zoom(ZOOM_FAC)

subplot(221)
hold on
draw_constraints(m, fids_list, v0);
draw_direction_field(m, DF2, [], scale_factor, cf_colors, face_alpha, edge_alpha, periods2_real)
%title({['MIQ p\in\mathbb{R}'], ['E = ', num2str(E2_real)]})
title({'MIQ real p', ['|AX-b|_2^2 = ', num2str(E2_real, '%.4f'), ...
                ', |p|_2^2 = ', num2str(norm(periods2_real)^2, '%.4f')]})
view([0,0,1]); axis off;
hold off
zoom(ZOOM_FAC)
% subplot(423)
% m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0);
% m.labelFaces(thetas2_real, '%.2f');
% view([0,0,1]); axis off;
% zoom(ZOOM_FAC)

subplot(224)
hold on
m.draw(ones(m.nV, 1), 'EdgeAlpha', 1.0, 'FaceAlpha', 0.0)
draw_local_frames(m, local_frames2)
m.draw(ones(m.nV, 1), 'FaceAlpha', 0.0, 'EdgeColor', period_tags2 > -1)
title('Local Frames'); hold off; view([0,0,1]); axis off
hold off

%% Plots without labels
figure(4); clf(4)
cf_colors = {'k', 'k', 'k', 'k'};
face_alpha = 1.0;
edge_alpha = 0.2;
scale_factor = m.avg_length;

subplot(223)
hold on
draw_constraints(m, fids_list, v0);
draw_direction_field(m, DF1, find(S1~=0), scale_factor, cf_colors, face_alpha, edge_alpha)
title({title1, ['|AX-b|_2^2 = ', num2str(E1, '%.4f'), ...
                ', |p|_2^2 = ', num2str(norm(periods1)^2, '%.4f')]})
view([0,0,1]); axis off;
hold off
zoom(ZOOM_FAC)

subplot(222)
hold on
draw_constraints(m, fids_list, v0);
draw_direction_field(m, DF2, [], scale_factor, cf_colors, face_alpha, edge_alpha)
%title({title2, ['E = ', num2str(E2)]})
title({title2, ['|AX-b|_2^2 = ', num2str(E2, '%.4f'), ...
                ', |p|_2^2 = ', num2str(norm(periods2)^2, '%.4f')]})
view([0,0,1]); axis off;
hold off
zoom(ZOOM_FAC)

subplot(221)
hold on
draw_constraints(m, fids_list, v0);
draw_direction_field(m, DF2, [], scale_factor, cf_colors, face_alpha, edge_alpha)
%title({['MIQ p\in\mathbb{R}'], ['E = ', num2str(E2_real)]})
title({'MIQ real p', ['|AX-b|_2^2 = ', num2str(E2_real, '%.4f'), ...
                ', |p|_2^2 = ', num2str(norm(periods2_real)^2, '%.4f')]})
view([0,0,1]); axis off;
hold off
zoom(ZOOM_FAC)

subplot(224)
hold on
draw_constraints(m, fids_list, v0);
S5 = find(d0'*x5 > 0.1);
for j = 1:size(S5,1)
    fid = S5(j, 1);
    H = plot3( m.V(fid,1), m.V(fid,2), m.V(fid,3), 'g.' );
    set( H, 'MarkerSize', 40 );
end
draw_direction_field(m, DF5, find(S1~=0), scale_factor, cf_colors, face_alpha, edge_alpha)

title(title5)
view([0,0,1]); axis off;
hold off
zoom(ZOOM_FAC)