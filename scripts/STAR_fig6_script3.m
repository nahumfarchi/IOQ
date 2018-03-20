p = find_data_folder();
FNAME = 'squares5_row.off';
fp = fullfile(p, FNAME);

N = 4;
%f0 = [1, 3];
%f0 = [1, 7];
f0 = [1, 19];
fc1 = f0(1);
fc2 = f0(2);
%v0 = normalize_rows([0, 1, 0; 1, 1, 0]);
%tc1 = pi+0.1;
%tc2 = pi-0.1;
tc1 = pi / 2 - 1e-5;
tc2 = pi / 4;
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

%% (4) min ||\eta||_1+||d0\theta+2\pi\eta||_2^2
%      s.t. \eta \in \mathbb{Z}
%           -\pi \leq \theta_i \leq \theta
%           \theta_i = \hat{\theta_i} \forall v_i \in V_c
d0_dual = zeros(m.niE, m.nF);
i = 1;
b = zeros(2*m.niE, 1);
for eid = 1:m.nE
    if m.isBoundaryEdge(eid)
        continue;
    end
    fi = m.EFAdj(eid, 1);
    fj = m.EFAdj(eid, 2);
    
    if fj == fc1
        b(i) = -tc1;
    elseif fj == fc2
        b(i) = -tc2;
    else
        d0_dual(i, fj) = 1;
    end
    
    if fi == fc1
        b(i) = tc1;
    elseif fi == fc2
        b(i) = tc2;
    else
        d0_dual(i, fi) = -1;
    end
    
    i = i + 1;
end

theta4 = zeros(m.nF, 1);
eta4 = zeros(m.niE, 1);
A0 = [d0_dual, zeros(m.niE, m.niE); zeros(m.niE, m.nF), -2*pi*eye(m.niE)];
A.times = @(x) A0*x;
A.trans = @(y) (y'*A0)';

opts.nonorth = true;
opts.tol = 1e-3;
opts.rho = 1e-3;
opts.weights = [zeros(m.nF, 1); ones(m.niE, 1)];

[x4, out] = yall1(A, b, opts);


%% Plots
ZOOM_FAC = 1.4;

figure(1); clf(1)
cf_colors = {'k', 'k', 'k', 'k'};
face_alpha = 0.0;
edge_alpha = 1.0;

subplot(425)
hold on
draw_constraints(m, f0, v0);
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
draw_constraints(m, f0, v0);
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
draw_constraints(m, f0, v0);
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


subplot(248)
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

%% Plots
figure(3); clf(3)
cf_colors = {'k', 'k', 'k', 'k'};
face_alpha = 0.0;
edge_alpha = 1.0;

subplot(223)
hold on
draw_constraints(m, f0, v0);
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
draw_constraints(m, f0, v0);
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
draw_constraints(m, f0, v0);
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

