% ========================================================================
% IOQ vs MIQ vs GO
% Run the three methods and create a simple plot that shows the energy,
% singularities, and timing.
% ========================================================================
FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = true;
SAVE = false;
mkdir('results')

%fp = '../../../data/torus_s0.off';
%fp = '../../../data/torus_fat_r2.off';
fp = '../../../data/ashish_nob/armadillo.off';
%fp = '../../../data/rounded_cube_keenan.off';
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;

% run IOQ
[alpha, beta, elapsed_ioq] = IOQ_highgenus_gpu(...
    V, F, ...
    'UseGPU', USE_GPU, ...
    'Iterations', 2000, ...
    'InvMethod', 'GPUBlockInv');

k = [alpha; beta];
res_ioq = TCODS(m, ...
    'k', k, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);
%[E_edges, Emiq] = per_edge_energy(res_tc);

% run MIQ
[res_miq, elapsed_miq] = ...
    nrosy_mex(fp, FACE0, GVEC, DEGREE);

% run GO
[res_go, Ego, elapsed_go] = GO(fp, 4);

R1 = res_ioq.ffield_vectors(1:nf, :);
R2 = res_miq.ffield_vectors(1:nf, :);
R3 = res_go.ffield_vectors(1:nf, :);
theta1 = res_ioq.ffield_angles;
theta2 = res_miq.ffield_angles;
theta3 = res_go.ffield_angles;

%% Plot        
title_ioq = {'IOQ', ...
    sprintf('Emiq = %.4g', res_ioq.miq_energy), ...
    sprintf('# sing = %d', res_ioq.n_vert_sing), ...
    sprintf('T = %.4g', elapsed_ioq)};
title_miq = {'MIQ', ...
    sprintf('Emiq = %.4g', res_miq.miq_energy), ...
    sprintf('# sing = %d', res_miq.n_vert_sing), ...
    sprintf('T = %.4g', elapsed_miq)};
title_go = {'GO', ...
    sprintf('Emiq = %.4g', res_go.miq_energy), ...
    sprintf('Ego = %.4g', Ego), ...
    sprintf('# sing = %d', res_go.n_vert_sing), ...
    sprintf('T = %.4g', elapsed_go)};
            
figure
opts = {'EdgeAlpha', 0.2, 'FaceAlpha', 1, 'PlotField', true};
subplot(131); res_ioq.draw(opts{:}); title(title_ioq)
subplot(132); res_miq.draw(opts{:}); title(title_miq)
subplot(133); res_go.draw(opts{:}); title(title_go)

if SAVE
    [~, meshname, ~] = fileparts(fp);
    filename = fullfile('results', meshname);
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end