%% ========================================================================
%  Time various Lx=b solvers
%  ========================================================================


data_folder = '../../../data/bunnies_large';

[~, n, ~] = fileparts(data_folder);
out_folder = fullfile('results', n);
if ~exist(out_folder, 'dir')
    mkdir(out_folder);
end
filepaths = get_filepaths(data_folder);
n_files = numel(filepaths);

SAVE = false;
REPS = 1;
SAME_SEED = true;
SEED = 112;
USE_GPU = false;
EPS = 0.5;
DROPTOL = 1e-4;
USE_GPU = false;
K = 1000;
BLOCK_SIZE = 1;
GD = gpuDevice();

TSetup = [];
TSolve = [];
TTotal = [];
ERR = [];
N_FACES = [];
names = {};
rng(SEED);

progressbar
for i = 1:1
    counter = 1;
    fp = filepaths{i};
    print_header(fp);

    m = Mesh(fp);
    nv = m.nV;
    d0 = get_exterior_derivatives(m);

    p = colamd(d0);
    clear pt;
    pt(p) = 1:length(p);
    tColamd = timeit(@() colamd(d0));
    d0p = d0(:, p);
    A = d0'*d0;
    Ap = d0p'*d0p;
    
    B = rand(nv, K);
    for j = 1:K
        B(:, j) = B(:, j) - mean(B(:, j));
    end
    
    % ---------------------------------------------------------------------
    % chol amd + pcg
    % ---------------------------------------------------------------------
    disp('ichol amd + pcg')
    tic
    Ri = ichol(Ap,struct('type','ict','droptol',DROPTOL,'michol','off'));
    tSetup = toc;

    tic
    
    %Ap = full(gpuArray(Ap));
    %B = gpuArray(B);
    %precond = full(gpuArray(Ri*Ri'));
    X = zeros(nv, K);
    disp('loop...')
    for j = 1:K
        if mod(j, 1) == 0
            disp(j);
        end
        c = B(p, j);
        %z = Ri' \ (Ri \ c);
        [z, flag] = pcg(Ap, c, 1e-6, 100, Ri, Ri');
        X(:, j) = z;
    end
    X = X(pt, :);
    
    %B = gather(B);
    %Ap = gather(Ap);
    
    wait(GD); tSolve = toc;
    tTotal = tSetup + tSolve;
    err = norm(A*X-B, 'fro') / norm(B, 'fro');
    
    TSetup(i, counter) = tSetup;
    TSolve(i, counter) = tSolve;
    TTotal(i, counter) = tTotal;
    ERR(i, counter) = err;
    counter = counter + 1;
    if i == 1
        names{end+1} = 'chol+amd+BS';
    end
    
    progressbar(i / n_files)
end
progressbar(1)