%% ========================================================================
%  Run JL IOQ on rounded cube with various eps and visualize the
%  singularities
%  ========================================================================

%% Setup
%fp = '../../../data/bunny.off';
%fp = '../../../data/round_cuber.off';
fp = '../../../data/rounded_cube_keenan.off';

FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0];

out_folder = fullfile('results');
mkdir(out_folder);
%if ~exist(out_folder, 'dir')
%    mkdir(out_folder);
%end
[~, filename, ~] = fileparts(fp);
filename = fullfile(out_folder, filename);

SAVE = true;
REPS = 30;
SAME_SEED = true;
SEED = 112;
USE_GPU = false;

m = Mesh(fp);
V = m.F; F = m.F; nv = m.nV; ne = m.nE;

results = {};
names = {};

%% ------------------------------------------------------------------------
%  True resistance
%% ------------------------------------------------------------------------
if SAME_SEED, rng(SEED); end
[alpha_p1, beta_p1, elapsed_total1] = IOQ_highgenus_gpu(V, F, ...
                        'highg_method', 'genus0', ...
                        'InvMethod', 'CholMexInv', ...
                        'bsx', false, ...
                        'Mesh', m, ...
                        'UseGPU', USE_GPU);
                    
k = [alpha_p1; beta_p1];
results{end+1} = TCODS(m, ...
    'k', k, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'gConstraintVec', GVEC);
names{end+1} = {'R', ['E=', num2str(results{end}.miq_energy)]};

%% ------------------------------------------------------------------------
%  Approx resistance
%% ------------------------------------------------------------------------
for eps = 0.3:0.1:1
    fprintf('eps = %g\n', eps);
    if SAME_SEED, rng(SEED); end
    [alpha_p2, beta_p2, elapsed_total2] = IOQ_highgenus_gpu(V, F, ...
        'InvMethod', 'ApproxResistance', ...
        'highg_method', 'genus0', ...
        'Iterations', 30, ...
        'UseGPU', USE_GPU, ...
        'Mesh', m, ...
        'JLEps', eps);
    
    k = [alpha_p2; beta_p2];
    results{end+1} = TCODS(m, ...
        'k', k, ...
        'f0', FACE0, ...
        'theta0', THETA0, ...
        'degree', DEGREE, ...
        'CreateFField', true, ...
        'gConstraintVec', GVEC);
    names{end+1} = {['Rtilde, eps=', num2str(eps)], ['E=',num2str(results{end}.miq_energy)]};
end

%% ------------------------------------------------------------------------
%  MIQ
%% ------------------------------------------------------------------------
%[theta, p, kmiq, R, local_frames, frame_diffs, elapsed, pFixed] = ...
%    NRoSy_mex_bin(fp, FACE0, GVEC, DEGREE);
%Emiq = E_MIQ(m, theta, frame_diffs, p, DEGREE);
[m, elapsed] = nrosy_mex(fp, FACE0, GVEC, DEGREE);
results{end+1} = m;
names{end+1} = {'MIQ', ['E=',num2str(results{end}.miq_energy)]};

%% ------------------------------------------------------------------------
%  Plot results
%% ------------------------------------------------------------------------
plotting_defaults(30, 30);
n_splots = length(results);
figure
for i = 1:n_splots 
    subplot(2, n_splots/2, i)
    results{i}.draw
    title(names{i})    
    view(3)
end
if SAVE
     %filename = fullfile(out_folder, 'ioq_cmp_timing_eps');
     print(filename, '-dpng', '-r300')
     saveas(gcf, filename, 'fig')
end

% %%
% plotting_defaults;
% K = 1;
% set(0,'DefaultAxesColorOrder',cbrewer('seq','YlOrRd',numel(names)-K));
% set(0,'DefaultFigureColormap',cbrewer('seq','YlOrRd',64));
% set2colors = cbrewer('qual', 'Set2', 8);
% 
% % energy
% figure(2); hold on
% [xx, I] = sort(N_FACES);
% Emean = mean(E, 3);
% Estd = std(E, 0, 3);
% 
% for i = 1:numel(names)-K
%     yy = Emean(I, i);
%     ee = Estd(I, i);
%     %plot(xx, yy, '-x')
%     errorbar(xx, yy, ee)
% end
% 
% yy = Emean(I, end);
% ee = Estd(I, end);
% %plot(xx, yy, '-x', 'Color', set2colors(1,:))
% errorbar(xx, yy, ee, 'Color', set2colors(1,:))
% 
% legend(names, 'interpreter', 'none', 'Location', 'northwest')
% title('IOQ final energy with R vs Rtilde')
% hold off
% 
% if SAVE
%     filename = fullfile(out_folder, 'IOQ_cmp_energy_eps');
%     print(filename, '-dpng', '-r300')
%     saveas(gcf, filename, 'fig')
% end
% 
% % timing
% figure(3); hold on
% [xx, I] = sort(N_FACES);
% Tmean = mean(T, 3);
% Tstd = std(T, 0, 3);
% 
% for i = 1:numel(names)-K
%     yy = Tmean(I, i);
%     ee = Tstd(I, i);
%     %plot(xx, yy, '-x')
%     errorbar(xx, yy, ee)
% end
% 
% yy = Tmean(I, end);
% ee = Tstd(I, end);
% %plot(xx, yy, '-x', 'Color', set2colors(1,:))
% errorbar(xx, yy, ee, 'Color', set2colors(1,:))
% 
% legend(names, 'interpreter', 'none', 'Location', 'northwest')
% title('IOQ timing with R vs Rtilde')
% hold off
% 
% if SAVE
%     filename = fullfile(out_folder, 'ioq_cmp_timing_eps');
%     print(filename, '-dpng', '-r300')
%     saveas(gcf, filename, 'fig')
% end
% 
% % save workspace
% if SAVE
%     filename = fullfile(out_folder, 'res');
%     save(filename, 'E', 'T', 'names');
% end
