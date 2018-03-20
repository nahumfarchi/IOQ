%%
%fp = '../../../data/bunny.off';
FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0];

%filepaths = get_filepaths('../../../data/bunnies');
data_folder = '../../../data/bunnies_large';
%data_folder = '../../../data/bunnies_small';
%data_folder = '../../../data/bunnies';
[~, n, ~] = fileparts(data_folder);
out_folder = fullfile('results', n);
if ~exist(out_folder, 'dir')
    mkdir(out_folder);
end
filepaths = get_filepaths(data_folder);
n_files = numel(filepaths);

SAVE = false;
REPS = 30;
SAME_SEED = false;
SEED = 112;

E = [];
T = [];
N_FACES = [];
names = {};
rng(SEED);
progressbar
prog_count = 0;
%%
%for i = 1:n_files
for i = 3:3
    
    fp = filepaths{i};
    %fp = '../../../data/bunny.off';
    print_header(fp);
    fprintf('%d // %d\n', i, n_files);
    m = Mesh(fp); V = m.V; F = m.F; ne = m.nE;

    for r = 1:REPS
    
        if i == 1 && r == 1
            names{end+1} = 'true_res';
        end
        j = 1;
        % True resistance
        print_header('True resistance');
        try
            if SAME_SEED, rng(SEED); end
            disp('Trying with GPU inv...')
            [alpha_p1, beta_p1, elapsed_total1] = IOQ_highgenus_gpu(V, F, ...
                'highg_method', 'genus0', ...
                'bsx', false, ...
                'Mesh', m);
        catch ME
            if SAME_SEED, rng(SEED); end
            try
                fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
                disp('Failed with GPU inv, trying with block inv...')
                [alpha_p1, beta_p1, elapsed_total1] = IOQ_highgenus_gpu(V, F, ...
                    'highg_method', 'genus0', ...
                    'InvMethod', 'GPUBlockInv', ...
                    'bsx', false, ...
                    'Mesh', m);
            catch ME
                try
                    fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
                    disp('Failed with block inv, trying with chol mex inv...')
                    [alpha_p1, beta_p1, elapsed_total1] = IOQ_highgenus_gpu(V, F, ...
                        'highg_method', 'genus0', ...
                        'InvMethod', 'CholMexInv', ...
                        'bsx', false, ...
                        'Mesh', m);
                catch ME
                    fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
                    disp('Failed with chol mex inv, trying without gpu iter...')
                    [alpha_p1, beta_p1, elapsed_total1] = IOQ_highgenus_gpu(V, F, ...
                        'highg_method', 'genus0', ...
                        'InvMethod', 'CholMexInv', ...
                        'bsx', false, ...
                        'Mesh', m, ...
                        'UseGPU', false);
                end
            end
            %k1 = [alpha_p1; beta_p1];
            %res1 = TCODS(m, 'k', k1, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
        end
        k1 = [alpha_p1; beta_p1];
        res1 = TCODS(m, 'k', k1, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
        E(i, j, r) = res1.miq_energy;
        T(i, j, r) = elapsed_total1;
        j = j + 1;

        % Approx resistance
        print_header('Approx resistance');
        for eps = 0.3:0.1:1
            if i == 1 && r == 1
                names{end+1} = ['approx_res_', num2str(eps)];
            end
            fprintf('eps = %g\n', eps);
            failed = false;
            try
                if SAME_SEED, rng(SEED); end
                disp('Trying with GPU...')
                [alpha_p2, beta_p2, elapsed_total2] = IOQ_highgenus_gpu(V, F, ...
                    'InvMethod', 'ApproxResistance', ...
                    'highg_method', 'genus0', ...
                    'Iterations', 30, ...
                    'UseGPU', true, ...
                    'Mesh', m, ...
                    'bsx', false, ...
                    'JLEps', eps, ...
                    'Colamd', true);
            catch ME
                try
                    if SAME_SEED, rng(SEED); end
                    fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
                    disp('Failed with GPU, trying without...')
                    [alpha_p2, beta_p2, elapsed_total2] = IOQ_highgenus_gpu(V, F, ...
                        'InvMethod', 'ApproxResistance', ...
                        'highg_method', 'genus0', ...
                        'Iterations', 30, ...
                        'UseGPU', false, ...
                        'Mesh', m, ...
                        'JLEps', eps, ...
                        'Colamd', true);
                catch ME
                    fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
                    fprintf('Failed with eps=%g\n', eps);
                    failed = true;
                end
            end
            if ~failed
                k2 = [alpha_p2; beta_p2];
                res2 = TCODS(m, 'k', k2, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
                E(i, j, r) = res2.miq_energy;
                T(i, j, r) = elapsed_total2;
                j = j + 1;
            else
                E(i, j, r) = nan;
                T(i, j, r) = nan;
                j = j + 1;
            end
        end


        % Approx resistance with sampling
    %     rng(112)
    %     [alpha_p3, beta_p3, elapsed_total3] = IOQ_highgenus_gpu(V, F, ...
    %             'InvMethod', 'ApproxResistance', ...
    %             'highg_method', 'genus0', ...
    %             'Iterations', 10, ...
    %             'UseGPU', true, ...
    %             'SampleResistance', true, ...
    %             'NSamples', 1000, ...
    %             'Mesh', m);
    %     k3 = [alpha_p3; beta_p3];
    %     res3 = TCODS(m, 'k', k3, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);

        % MIQ
        if i == 1 && r == 1
            names{end+1} = 'MIQ';   
        end
        tic
        [theta, p, kmiq, R, local_frames, frame_diffs, elapsed, pFixed] = ...
            NRoSy_mex_bin(fp, FACE0, GVEC, DEGREE);
        elapsed_total4 = toc;
        Emiq = E_MIQ(m, theta, frame_diffs, p, DEGREE);
        E(i, j, r) = Emiq;
        T(i, j, r) = elapsed_total4;
        j = j + 1;

        %fprintf('res1 E = %g\nres2 E = %g\nres3 E = %g\nMIQ E = %g\n', res1.miq_energy, res2.miq_energy, res3.miq_energy, Emiq);
        res1.miq_energy
        res2.miq_energy
        Emiq

        %figure(1); 
        %subplot(121); res1.draw; title({'IOQ with true resistance', sprintf('E = %g', res1.miq_energy), sprintf('T = %g', elapsed_total1)})
        %subplot(122); res2.draw; title({'IOQ with approx resistance', sprintf('E = %g', res2.miq_energy), sprintf('T = %g', elapsed_total2)})
        prog_count = prog_count+1;
        progressbar(prog_count / (n_files*REPS))
    end
    
    
    
    N_FACES(end+1) = m.nF;
end

%%
plotting_defaults(30, 30);
K = 1;
set(0,'DefaultAxesColorOrder',cbrewer('seq','YlOrRd',numel(names)-K));
set(0,'DefaultFigureColormap',cbrewer('seq','YlOrRd',64));
set2colors = cbrewer('qual', 'Set2', 8);

%out_folder = fullfile('results', 'bunnies_small_ichol_amd_pcg');
%mkdir(out_folder);
%SAVE = true;

% energy
figure(2); hold on
[xx, I] = sort(N_FACES);
Emean = mean(E, 3);
Estd = std(E, 0, 3);

for i = 1:numel(names)-K
    yy = Emean(I, i);
    ee = Estd(I, i);
    %plot(xx, yy, '-x')
    errorbar(xx, yy, ee)
end

yy = Emean(I, end);
ee = Estd(I, end);
%plot(xx, yy, '-x', 'Color', set2colors(1,:))
errorbar(xx, yy, ee, 'Color', set2colors(1,:))

legend(names, 'interpreter', 'none', 'Location', 'northwest')
title('IOQ final energy with R vs Rtilde')
hold off

if SAVE
    filename = fullfile(out_folder, 'IOQ_cmp_energy_eps');
    print(filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

% timing
figure(3); hold on
[xx, I] = sort(N_FACES);
Tmean = mean(T, 3);
Tstd = std(T, 0, 3);

for i = 1:numel(names)-K
    yy = Tmean(I, i);
    ee = Tstd(I, i);
    %plot(xx, yy, '-x')
    errorbar(xx, yy, ee)
end

yy = Tmean(I, end);
ee = Tstd(I, end);
%plot(xx, yy, '-x', 'Color', set2colors(1,:))
errorbar(xx, yy, ee, 'Color', set2colors(1,:))

legend(names, 'interpreter', 'none', 'Location', 'northwest')
title('IOQ timing with R vs Rtilde')
hold off

if SAVE
    filename = fullfile(out_folder, 'ioq_cmp_timing_eps');
    print(filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

% save workspace
if SAVE
    filename = fullfile(out_folder, 'res');
    save(filename, 'E', 'T', 'names', 'N_FACES');
end

%%
% if ~exist('results', 'dir')
%     mkdir('results')
% end
% 
% figure(2);
% plot(N_FACES, E1, '--X', N_FACES, E2, '--X', N_FACES, E4, '--X');
% legend('R', 'Rtilde', 'MIQ')
% title('IOQ final energy with R vs Rtilde')
% if SAVE
%     print('results/IOQ_cmp_energy', '-dpng', '-r300')
% end
% 
% figure(3);
% plot(N_FACES, T1, '--X', N_FACES, T2, '--X', N_FACES, T4, '--X');
% legend('R', 'Rtilde', 'MIQ')
% title('IOQ elapsed time with R vs Rtilde')
% if SAVE
%     print('results/IOQ_cmp_timing', '-dpng', '-r300')
% end