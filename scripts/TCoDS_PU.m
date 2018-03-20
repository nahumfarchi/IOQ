%% Settings
%FNAME = 'sphere_s0.off';
%FNAME = 'sphere_s0.off';
FNAME = 'sphere_s2.off';
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

%[local_frames2, r2] = create_local_frames(m);


%% (4) TCoDS with libigl singularities (with a minus sign)
title4 = '(4) TCoDS with libigl singularities';
print_header(title4)
inds = find(abs(S1) > 1e-10);
k4 = [inds, -N*S1(inds)];
connection4 = create_trivial_connection(m, k4, N, VERBOSE);
E4 = norm(connection4)^2;
DF4 = create_direction_field(m, connection4, f0, v0, N);

y4 = r1+(pi/2)*periods1+d1'*thetas1;
lambda = d0 \ (2*y4);
norm(d0*lambda-2*y4)


%% Plots
figure(1); clf(1)
cf_colors = {'k', 'k', 'k', 'k'};

subplot(221)
draw_direction_field(m, DF1, [], scale, cf_colors)
subplot(222)
draw_direction_field(m, DF4, [], scale, cf_colors)