
% Compare the following: 
%
% MIQ --> (p,theta) --> k --> TCoDS --> x* --> (p_n, theta_n) -->
% check_norm(p - p_n, theta - theta_n)
%
% MIQ --> (p,theta) --> k --> change of vars to TCoDS --> x --> (p_n,
% theta-N) --> check_norm(p - p_n, theta - theta_n)
%
% random k --> TCoDS --> x* --> (p,theta) --> k_n --> check_norm(k - k_n)

%% Settings
%FNAME = 'sphere_s0.off';
%FNAME = 'sphere_s1.off';
%FNAME = 'ellipsoid_s4r.off';
%FNAME = 'round_cuber.off';
%FNAME = 'sphere_s0.off';
%FNAME = 'torus_fat_r2.off';
FNAME = 'sphere_s0.off';
VERBOSE = true;
EPS = 1e-10;

% Find mesh file
p = find_data_folder();
fp = fullfile(p, FNAME);

% Load mesh
m = Mesh();
m.loadTM(fp);
scale = 2 * m.avg_length;

[d0, d1] = get_exterior_derivatives(m);
%d0 = -d0;
%d1 = -d1;
N = 4;          % Direction field degree (i.e., 4 for cross field)
f0 = [1];       % Starting face
%v0 = [0, 0, 1]; % Starting direction
v0 = m.V(m.F(f0, 2), :) - m.V(m.F(f0, 1), :);


%% libigl MIQ (greedy rounding)
print_header('libigl MIQ')
res1 = nrosy_wrapper(fp, f0, v0, N);
theta0 = [res1.local_frames(1,:); res1.local_frames(m.nF+2,:)] * v0';
theta0 = atan2(theta0(2), theta0(1));

%% MIQ --> TCoDS
print_header('MIQ --> TCoDS')
Kg = res1.Kg;
local_frames = res1.local_frames;
r = res1.frame_diffs;
[my_local_frames, my_r] = create_local_frames(m);
my_r = wrapToPi(my_r);
p = res1.p;
theta = res1.theta;

[res2.x, res2.k] = MIQ_to_TCoDS(theta, p, r, Kg, d0, d1);
res2.nrosy = create_direction_field(m, res2.x, f0, theta0, N, local_frames, r);
res2.E = norm(res2.x)^2;

check_norm('res2.E', 'res1.E');
check_norm('d1*res2.x', '0');
%norm(res2.E - res1.E)
%norm(d1*res2.x)


%% TCoDS --> MIQ
print_header('TCoDS --> MIQ')
[res3.theta, res3.p] = TCoDS_to_MIQ(res2.x, res2.k, r, theta0, Kg, d0, d1, res1.period_is_fixed);
res3.E = E_MIQ(m, res3.theta, r, res3.p, N);

check_norm('res3.p', 'res1.p');
check_norm('res3.theta', 'res1.theta');
check_norm('res3.E', 'res1.E');
%norm(res3.p - res1.p)
%norm(res3.theta - res1.theta)


%% TCoDS with libigl singularities
print_header('TCoDS with libigl singularities')
%I = [res1.S_inds, -res1.S(res1.S_inds)];
inds = find(abs(res2.k) > EPS);
I = [inds, res2.k(inds)/N];
theta0 = [res1.local_frames(1,:); res1.local_frames(m.nF+1,:)] * v0';
theta0 = atan2(theta0(2), theta0(1));
[res4.x, res4.E] = create_trivial_connection(m, I, N, VERBOSE);

check_norm('res4.x', 'res2.x');
check_norm('d1*res4.x', '0');
check_norm('res4.E', 'res1.E');
%assert(norm(res4.x - res2.x) < EPS)
%norm(d1*res4.x)


%% opt TCoDS --> MIQ
print_header('opt TCoDS --> MIQ')
[res5.theta, res5.p] = TCoDS_to_MIQ(res4.x, res2.k, r, theta0, Kg, d0, d1, res1.period_is_fixed);
res5.E = E_MIQ(m, res5.theta, r, res5.p, N);

check_norm('res5.p', 'res1.p');
check_norm('res5.theta', 'res1.theta');
check_norm('res5.E', 'res1.E');


%% Invent a random k
print_header('Randon k')
xi = 2 - 2*m.genus;
res6.k = zeros(m.nV, 1);
for i = 1:N*xi
    ind = floor(rand*m.nV)+1;
    while res6.k(ind) > 0
        ind = floor(rand*m.nV)+1;
    end
    res6.k(ind) = 1;
end

inds = find(abs(res6.k) > EPS);
I = [inds, res6.k(inds)/N];
theta0 = [res1.local_frames(1,:); res1.local_frames(m.nF+1,:)] * v0';
theta0 = atan2(theta0(2), theta0(1));
[res6.x, res6.E] = create_trivial_connection(m, I, N, VERBOSE);

[res6.theta, res6.p] = TCoDS_to_MIQ(res6.x, res6.k, r, theta0, Kg, d0, d1, res1.period_is_fixed);
[res6.x_n, res6.k_n] = MIQ_to_TCoDS(res6.theta, res6.p, r, Kg, d0, d1);

check_norm('d1*res6.x', '0');
check_norm('d1*res6.x_n', '0');
check_norm('res6.k', 'res6.k_n');
