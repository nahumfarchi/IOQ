%% Settings
FNAME = 'sphere_s1.off';
%FNAME = 'bumpy.off';
%FNAME = 'bunny.off';
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
f0 = [1, 100, 200, 300, 400, 500];
v0 = [m.V(m.F(1, 2), :) - m.V(m.F(1, 1), :);
      m.V(m.F(100, 2), :) - m.V(m.F(100, 1), :);
      m.V(m.F(200, 2), :) - m.V(m.F(200, 1), :);
      m.V(m.F(300, 2), :) - m.V(m.F(300, 1), :);
      m.V(m.F(400, 2), :) - m.V(m.F(400, 1), :);
      m.V(m.F(500, 2), :) - m.V(m.F(500, 1), :)];
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

[local_frames2, r2] = create_local_frames(m);


%% (2) MIQ --> TCoDS change of variables
y2 = r1 + (pi/2)*periods1 + d1'*thetas1;
norm(d1' * ( inv(d1*d1') * d1 * y2))
norm(d1' * ( (d1*d1') \ (d1*y2)))
norm(d1' * ( (d1*d1') \ d1 * y2))