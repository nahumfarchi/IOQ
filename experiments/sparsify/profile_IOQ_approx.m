%%
fp = '../../../data/bunnies/bunny_26k_faces.off';
m = Mesh(fp);
ne = m.nE;

%% true resistance

profile clear
profile on

try
    rng(112)
    disp('Trying with GPU inv...')
    [alpha_p1, beta_p1, elapsed_total1] = IOQ_highgenus_gpu_BAKK(V, F, ...
        'highg_method', 'genus0', ...
        'bsx', true, ...
        'Mesh', m);
    %k1 = [alpha_p1; beta_p1];
    %res1 = TCODS(m, 'k', k1, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
catch
    rng(112)
    disp('Failed with GPU inv, trying with block inv...')
    [alpha_p1, beta_p1, elapsed_total1] = IOQ_highgenus_gpu_BAKK(V, F, ...
        'highg_method', 'genus0', ...
        'InvMethod', 'GPUBlockInv', ...
        'bsx', false, ...
        'Mesh', m);
    %k1 = [alpha_p1; beta_p1];
    %res1 = TCODS(m, 'k', k1, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
end

profile viewer

%% approx resistance

profile clear
profile on

try
    rng(112)
    disp('Trying with GPU...')
    [alpha_p2, beta_p2, elapsed_total2] = IOQ_highgenus_gpu_BAKK(V, F, ...
        'InvMethod', 'ApproxResistance', ...
        'highg_method', 'genus0', ...
        'Iterations', 30, ...
        'UseGPU', true, ...
        'Mesh', m);
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

profile viewer