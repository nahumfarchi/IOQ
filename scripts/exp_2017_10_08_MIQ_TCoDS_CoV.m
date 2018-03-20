%% Setup
fp = '../data/bunny.off';
%fp = '../data/torus_fat_r2.off';
%vert_sing = [...
%    1, -0.5; ...
%    100, -0.5; ...
%    200, 1.5; ...
%    300, 0.5; ...
%    400, 1;];
f0 = 1;
theta0 = 0;
degree = 4;
use_cotan = false;
constrained_faces = [1];
constraint_vects = [1, 0, 0];

m = Mesh();
m.loadTM(fp);

disp('---')


%%
print_header('mex MIQ')
tic
res_miq = nrosy_mex(fp, constrained_faces, constraint_vects, degree);
toc

%theta0 = [miq.local_frames(1,:); miq.local_frames(m.nF+2,:)] * constraint_vects(1,:)';
%theta0 = atan2(theta0(2), theta0(1));

%%
print_header('mex tcods with miq singularities')
tic
[ffield, E, x] = tcods_mex(fp, res_miq.S(:,1), degree*res_miq.S(:,2), [], [], f0, theta0, degree, use_cotan);
%[ffield, theta]  = connection_to_nrosy(m, x, f0, theta0, degree);
toc

res_tcods_crane.ffield = ffield;
res_tcods_crane.degree = degree;
res_tcods_crane.S = res_miq.S;


%%
print_header('matlab tcods with miq singularities')
tic
tcods_mine = TCODS(m, res_miq.S, f0, theta0, degree, false);
toc

%%
print_header('MIQ(theta, p) --> TCoDS(x, k)')
[d0, d1] = get_exterior_derivatives(m);
Kg = get_gaussian_curvature(m);
%Kg = res1.Kg;
local_frames = res_miq.local_frames;
r = res_miq.frame_diffs;
%[my_local_frames, my_r] = create_local_frames(m);
%my_r = wrapToPi(my_r);
p = res_miq.p;
theta = res_miq.theta;

[res_miq_to_tcods.x, res_miq_to_tcods.k] = ...
    MIQ_to_TCoDS(theta, p, r, Kg, d0, d1, degree);
%res_miq_to_tcods.ffield = create_direction_field(m, res_miq_to_tcods.x, f0, theta0, degree, local_frames, r);
[res_miq_to_tcods.ffield, res_miq_to_tcods.theta] = ...
    connection_to_nrosy(m, res_miq_to_tcods.x, f0, theta0, degree, 'LocalFrames', local_frames, 'FrameDiffs', r);
res_miq_to_tcods.E = norm(res_miq_to_tcods.x)^2;
res_miq_to_tcods.degree = degree;

inds = find(abs(res_miq_to_tcods.k) > 1e-10);
res_miq_to_tcods.S = [inds, res_miq_to_tcods.k(inds)];

check_norm('res_miq_to_tcods.E', 'res_miq.E');
check_norm('d1*res_miq_to_tcods.x', '0');

%%
print_header('TCoDS(x, k) --> MIQ(theta, p)')
[res_tcods_to_miq.theta, res_tcods_to_miq.p] = ...
    TCoDS_to_MIQ(res_miq_to_tcods.x, res_miq_to_tcods.k, r, theta0, Kg, d0, d1, res_miq.pFixed);
res_tcods_to_miq.E = ...
    E_MIQ(m, res_tcods_to_miq.theta, r, res_tcods_to_miq.p, degree);
res_tcods_to_miq.ffield = angles_to_ffield(res_tcods_to_miq.theta, local_frames, degree);
res_tcods_to_miq.degree = degree;


%% Plot
figure()
subplot(231)
MeshVis.plot(m, 'nrosy', tcods_mine);
title(['My tcods with miq singularities, E = ', num2str(tcods_mine.E)])

subplot(232)
MeshVis.plot(m, 'nrosy', res_tcods_crane);
title(['Crane''s tcods with miq singularities, E = ', num2str(E)])

subplot(233)
MeshVis.plot(m, 'nrosy', res_miq);
title(['MIQ, E = ', num2str(res_miq.E)])

subplot(234)
MeshVis.plot(m, 'nrosy', res_miq_to_tcods)
title(['miq(theta, p) --> tcods(x, k), E = ', num2str(res_miq_to_tcods.E)])

subplot(235)
MeshVis.plot(m, 'nrosy', res_tcods_to_miq)
title(['tcods(x, k) --> miq(theta, p), E = ', num2str(res_tcods_to_miq.E)])


