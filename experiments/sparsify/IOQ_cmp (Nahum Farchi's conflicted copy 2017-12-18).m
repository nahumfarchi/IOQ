%fp = '../../../data/bunny.off';
FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0];

filepaths = get_filepaths('../../../data/bunnies');
n_files = numel(filepaths);

SAVE = false;

E = [];
T = [];
N_FACES = [];
for i = 1:n_files
    fp = filepaths{i};
    %fp = '../../../data/bunny.off';
    print_header(fp);
    fprintf('%d // %d\n', i, n_files);
    m = Mesh(fp); V = m.V; F = m.F; ne = m.nE;

    j = 1;
    % True resistance
    try
        rng(112)
        disp('Trying with GPU inv...')
        [alpha_p1, beta_p1, elapsed_total1] = IOQ_highgenus_gpu_BAKK(V, F, ...
            'highg_method', 'genus0', ...
            'bsx', false, ...
            'Mesh', m);
    catch
        rng(112)
        try
            disp('Failed with GPU inv, trying with block inv...')
            [alpha_p1, beta_p1, elapsed_total1] = IOQ_highgenus_gpu_BAKK(V, F, ...
                'highg_method', 'genus0', ...
                'InvMethod', 'GPUBlockInv', ...
                'bsx', false, ...
                'Mesh', m);
        catch
            try
                disp('Failed with block inv, trying with chol mex inv...')
                [alpha_p1, beta_p1, elapsed_total1] = IOQ_highgenus_gpu_BAKK(V, F, ...
                    'highg_method', 'genus0', ...
                    'InvMethod', 'CholMexInv', ...
                    'bsx', false, ...
                    'Mesh', m);
            catch
                disp('Failed with chol mex inv, trying without gpu iter...')
                [alpha_p1, beta_p1, elapsed_total1] = IOQ_highgenus_gpu_BAKK(V, F, ...
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
    E(i, j) = res1.miq_energy;
    T(i, j) = elapsed_total1;
    j = j + 1;
    
    % Approx resistance
    for eps = 0.1:0.1:1
        try
            rng(112)
            disp('Trying with GPU...')
            [alpha_p2, beta_p2, elapsed_total2] = IOQ_highgenus_gpu_BAKK(V, F, ...
                'InvMethod', 'ApproxResistance', ...
                'highg_method', 'genus0', ...
                'Iterations', 30, ...
                'UseGPU', true, ...
                'Mesh', m, ...
                'bsx', false);
        catch
            rng(112)
            disp('Failed with GPU, trying without...')
            [alpha_p2, beta_p2, elapsed_total2] = IOQ_highgenus_gpu_BAKK(V, F, ...
                'InvMethod', 'ApproxResistance', ...
                'highg_method', 'genus0', ...
                'Iterations', 30, ...
                'UseGPU', false, ...
                'Mesh', m);
        end
        k2 = [alpha_p2; beta_p2];
        res2 = TCODS(m, 'k', k2, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
        E(i, j) = res2.miq_energy;
        T(i, j) = elapsed_total2;
        j = j + 1;
    end
    
    
    % Approx resistance with sampling
%     rng(112)
%     [alpha_p3, beta_p3, elapsed_total3] = IOQ_highgenus_gpu_BAKK(V, F, ...
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
    tic
    [theta, p, kmiq, R, local_frames, frame_diffs, elapsed, pFixed] = ...
        NRoSy_mex_bin(fp, FACE0, GVEC, DEGREE);
    elapsed_total4 = toc;
    Emiq = E_MIQ(m, theta, frame_diffs, p, DEGREE);
    E(i, j) = Emiq;
    T(i, j) = elapsed_total4;
    j = j + 1;

    %fprintf('res1 E = %g\nres2 E = %g\nres3 E = %g\nMIQ E = %g\n', res1.miq_energy, res2.miq_energy, res3.miq_energy, Emiq);
    res1.miq_energy
    res2.miq_energy
    Emiq

    figure(1); 
    subplot(121); res1.draw; title({'IOQ with true resistance', sprintf('E = %g', res1.miq_energy), sprintf('T = %g', elapsed_total1)})
    subplot(122); res2.draw; title({'IOQ with approx resistance', sprintf('E = %g', res2.miq_energy), sprintf('T = %g', elapsed_total2)})
    
    N_FACES(end+1) = m.nF;
end

[N_FACES, I] = sort(N_FACES);
E1 = E1(I);
T1 = T1(I);
E2 = E2(I);
T2 = T2(I);
%E3 = E3(I);
%T3 = T3(I);
E4 = E4(I);
T4 = T4(I);

if ~exist('results', 'dir')
    mkdir('results')
end

figure(2);
plot(N_FACES, E1, '--X', N_FACES, E2, '--X', N_FACES, E4, '--X');
legend('R', 'Rtilde', 'MIQ')
title('IOQ final energy with R vs Rtilde')
if SAVE
    print('results/IOQ_cmp_energy', '-dpng', '-r300')
end

figure(3);
plot(N_FACES, T1, '--X', N_FACES, T2, '--X', N_FACES, T4, '--X');
legend('R', 'Rtilde', 'MIQ')
title('IOQ elapsed time with R vs Rtilde')
if SAVE
    print('results/IOQ_cmp_timing', '-dpng', '-r300')
end