%% Settings
%FNAME = 'bumpy.off';
FNAME = 'sphere_s1.off';
OUT_FOLDER = '21.5.17';
VERBOSE = true;

mkdir(OUT_FOLDER);

% Find and load mesh file
DATA_FOLDER = find_data_folder();
fp = fullfile(DATA_FOLDER, FNAME);

m = Mesh();
m.loadTM(fp);

scale = 4 * m.avg_length;

N = 1;
f0 = [1]; %, 100];
v0 = [m.V(m.F(1, 2), :) - m.V(m.F(1, 1), :)];
%      m.V(m.F(100, 2), :) - m.V(m.F(100, 1), :)];
%v0 = v0 / norm(v0);
v0 = normalize_rows(v0);

[d0, d1] = get_exterior_derivatives(m);


%% (1) libigl MIQ greedy rounding
title1 = '(1) libigl MIQ greedy rounding';
print_header(title1)
[R1, S1, periods1, period_is_fixed1, thetas1, local_frames1, r1, EV1, EF1, FE1] = nrosy_wrapper(fp, f0, v0, N);
E1 = E_MIQ(m, thetas1, r1, periods1, N);
DF1 = angles_to_DF(thetas1, local_frames1, N);
disp(['E1 = ', num2str(E1)])


%% (2) MIQ system with libigl periods
title2 = '(2) MIQ system with libigl periods';
print_header(title2)
[A2, b2, thetas2, periods2, theta_tags2, period_tags2] = create_MIQ_system(m, N, f0, v0, local_frames1, r1, periods1, ones(size(periods1))); 
[thetas2, periods2, E2] = solve_MIQ(A2, b2, thetas2, periods2, theta_tags2, period_tags2);
DF2 = angles_to_DF(thetas2, local_frames1, N);
disp(['E2 = ', num2str(E2)])


%% (3) TCoDS with k=d_0^T p, where p are the periods given by libigl greedy rounding.
title3 = '(3) TCoDS (k=d_0^T p)';
print_header(title3)
%y3 = d1' * thetas2 + r1 + 2*pi/N * periods2;
y3 = d1'*thetas2 + r1 + pi/2 * periods2;
disp(['y3 = ', num2str(norm(y3)^2)])
%S3 = (get_gaussian_curvature(m) + d0'*r1) / (2*pi) + d0'*periods2/4;
k3 = d0' * periods2;
k3 = [(1:length(k3))', k3];
%inds = find(abs(S3) > 1e-10);
%k3 = [inds, -S3(inds)];
connection3 = create_trivial_connection(m, k3, N, VERBOSE);
E3 = norm(connection3)^2;
DF3 = create_direction_field(m, connection3, f0, v0, N);



%% (4) TCoDS with libigl singularities (with a minus sign)
title4 = '(4) TCoDS with libigl singularities';
print_header(title4)
inds = find(abs(S1) > 1e-10);
k4 = [inds, -N*S1(inds)];
connection4 = create_trivial_connection(m, k4, N, VERBOSE);
E4 = norm(connection4)^2;
DF4 = create_direction_field(m, connection4, f0, v0, N);


%% (5) MIQ system direct rounding
title5 = '(5) MIQ system direct rounding';
print_header(title5)
[A5, b5, thetas5, periods5, theta_tags5, period_tags5] = create_MIQ_system(m, N, f0, v0, local_frames1, r1);
[thetas5, periods5, E5] = solve_MIQ(A5, b5, thetas5, periods5, theta_tags5, period_tags5);
disp(['E5 = ', num2str(E5)])
DF2 = angles_to_DF(thetas5, local_frames1, N);

%% (6) TCoDSk with some arbitrary set of singularities
title6 = '(6) TCoDS with some arbitrary set of singularities';
print_header(title6)
k6 = [300, 1; 10 1];
connection6 = create_trivial_connection(m, k6, N, VERBOSE);
E6 = norm(connection6)^2;
DF6 = create_direction_field(m, connection6, f0, v0, N);


%% Plots
figure(1)

subplot(221)
cf_colors = {'k', 'k', 'k', 'k'};
draw_direction_field(m, DF1, k4, scale/4, cf_colors)
title({title1, num2str(E1)})

subplot(222)
draw_direction_field(m, DF2, [], scale/4, cf_colors)
title({title2, num2str(E2)})

subplot(223)
draw_direction_field(m, DF3, [], scale/4, cf_colors)
title({title3, num2str(E3)})

subplot(224)
draw_direction_field(m, DF4, [], scale/4, cf_colors)
title({title4, num2str(E4)})

figure(2);
%ha(1) = subplot(131);
ha(1) = gca;
draw_direction_field(m, DF1, k4, 10*scale, cf_colors)
%title({'MIQ Greedy Rounding', num2str(E1)})
%title(['E=', num2str(E1)])

%ha(3) = subplot(132);
figure(3)
ha(2) = gca;
draw_direction_field(m, DF6, k6, scale/4, cf_colors)
%title({'TCoDS With Arbitrary Singularities', num2str(E6)})
title(['E=', num2str(E6)])

%ha(2) = subplot(133);
figure(4)
ha(3) = gca;
draw_direction_field(m, DF4, k4, scale/2, cf_colors, 1.0, 0.0)
%title({'TCoDS + MIQ Singularities', num2str(E4)})
%title(['E=', num2str(E4)])

hlink = linkprop(ha, {'CameraPosition', 'CameraUpVector'});