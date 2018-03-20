%% Setup
fp = '../data/torus_fat_r2.off';
%vert_sing = [...
%    1, -0.5; ...
%    100, -0.5; ...
%    200, 1.5; ...
%    300, 0.5; ...
%    400, 1;];
f0 = 1;
theta0 = pi / 2;
degree = 1;
use_cotan = false;
constrained_faces = [1];
constraint_vects = [1, 0, 0];

m = Mesh();
m.loadTM(fp);

disp('---')

% f0 = [1];
% v0 = [1,0,0];
% %NRoSy_mex_bin(fp, f0, v0, degree);
% res_miq = nrosy_mex(fp, constrained_faces, constraint_vects, degree);

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
[ffield, E, x] = tcods_mex(fp, res_miq.S(:,1), degree*res_miq.S(:,2), [1, 2], [0, 0], f0, theta0, degree, use_cotan);
%[ffield, theta]  = connection_to_nrosy(m, x, f0, theta0, degree);
toc

tcods0.ffield = ffield;
tcods0.degree = degree;
tcods0.S = res_miq.S;
tcods0.E = E;

%%
print_header('mex tcods with miq singularities')
tic
[ffield, E, x] = tcods_mex(fp, res_miq.S(:,1), degree*res_miq.S(:,2), [1, 2], [2, 0], f0, theta0, degree, use_cotan);
%[ffield, theta]  = connection_to_nrosy(m, x, f0, theta0, degree);
toc

tcods2pi.ffield = ffield;
tcods2pi.degree = degree;
tcods2pi.S = res_miq.S;
tcods2pi.E = E;

%%
print_header('mex tcods with miq singularities')
tic
[ffield, E, x] = tcods_mex(fp, res_miq.S(:,1), degree*res_miq.S(:,2), [1, 2], [4, 0], f0, theta0, degree, use_cotan);
%[ffield, theta]  = connection_to_nrosy(m, x, f0, theta0, degree);
toc

tcods4pi.ffield = ffield;
tcods4pi.degree = degree;
tcods4pi.S = res_miq.S;
tcods4pi.E = E;

%%
print_header('mex tcods with miq singularities')
tic
[ffield, E, x] = tcods_mex(fp, res_miq.S(:,1), degree*res_miq.S(:,2), [1, 2], [6, 0], f0, theta0, degree, use_cotan);
%[ffield, theta]  = connection_to_nrosy(m, x, f0, theta0, degree);
toc

tcods6pi.ffield = ffield;
tcods6pi.degree = degree;
tcods6pi.S = res_miq.S;
tcods6pi.E = E;


%%
print_header('matlab tcods with miq singularities')
tic
tcods_mine = TCODS(m, res_miq.S, f0, theta0, degree, false);
toc

%% Plot
figure()

subplot(221)
MeshVis.plot(m, 'nrosy', tcods0, 'Scale', 4)
title(['Generator holonomy: 0, E: ', num2str(tcods0.E)])
view(0, 0)

subplot(222)
MeshVis.plot(m, 'nrosy', tcods2pi, 'Scale', 4)
title(['Generator holonomy: 2pi, E: ', num2str(tcods2pi.E)])
view(0, 0)

subplot(223)
MeshVis.plot(m, 'nrosy', tcods4pi, 'Scale', 4)
title(['Generator holonomy: 4pi, E: ', num2str(tcods4pi.E)])
view(0, 0)

subplot(224)
MeshVis.plot(m, 'nrosy', tcods6pi, 'Scale', 4)
title(['Generator holonomy: 6pi, E: ', num2str(tcods6pi.E)])
view(0, 0)

%% 
figure()
MeshVis.plot(m, 'nrosy', tcods6pi, 'Scale', 4)
view(0, 0)