FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = false;
fp = '../data/round_cuber.off';
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;
CONSTRAINTS = [1, 0; 50, pi; 100, 0.1];
PLOT_PROPS = {'FaceColor', 'w', 'PlotField', true, 'Constraints', CONSTRAINTS};

%% run IOQ
[alpha, beta, x, stats, out] = IOQ(...
    V, F, ...
    'GPU', USE_GPU, ...
    'Iterations', 2000, ...
    'Constraints', CONSTRAINTS);
gamma = out.gamma;
disp(norm(x)^2)

%% run TCODS
k = [alpha; beta];

res_tc = TCODS(m, ...
    'k', k, ...
    'f0', FACE0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'Constraints', CONSTRAINTS, ...
    'Gamma', gamma);
%[E_edges, Emiq] = per_edge_energy(res_tc);

% TODO check with gamma
theta = res_tc.ffield_angles;

for i = 2:size(CONSTRAINTS, 1)
    fid = CONSTRAINTS(i, 1);
    t = CONSTRAINTS(i, 2);
    assert(norm(theta(fid) - (t + pi/2 * gamma(i-1))) < 1e-10)
    %fprintf('theta(%d) = %g\nt + pi/2*gamma(%d)=%g\n', fid, theta(fid), i-1, t+pi/2*gamma(i-1));
end

disp(res_tc.miq_energy)

%%
cam = [];
% cam.pba = [187.409091 342.300000 187.409091 ];
% cam.dar = [1 1 1 ];
% cam.cva = [2.181558 ];
% cam.cuv = [-0.046484 0.797630 -0.601353 ];
% cam.ct = [0.000080 0.000084 -0.000167 ];
% cam.cp = [2.209290 -13.807737 -18.485526 ];

title1 = sprintf('ioq, E = %g', out.m.miq_energy);
title2 = sprintf('tc,  E = %g', res_tc.miq_energy);

figure; out.m.draw(PLOT_PROPS{:}); title(title1)
set_camera(gca, cam);
figure; res_tc.draw(PLOT_PROPS{:}); title(title2)
set_camera(gca, cam);

hold on
%res_tc.drawLabels()
hold off