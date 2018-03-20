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
K = 1;
BLOCK_SIZE = 1;

TSetup = [];
TSolve = [];
TTotal = [];
ERR = [];
N_FACES = [];
names = {};
rng(SEED);
for i = 1:n_files
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
    
    %k = round(24*log(nv)/EPS^2);
    
    B = rand(nv, K);
    for j = 1:K
        B(:, j) = B(:, j) - mean(B(:, j));
    end

    % ---------------------------------------------------------------------
    % chol amd + BS
    % ---------------------------------------------------------------------
    disp('chol amd + BS')
    tic
    [R, ~, S] = chol(A);
    tSetup = toc;

    tic
    %%C = S'*B;
    X = S * (R \ (R' \ (S' * B)));
    %X = S * Z;
    tSolve = toc;
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
    
    %% --------------------------------------------------------------------
    % chol amd + parfor
    % ---------------------------------------------------------------------
%     disp('chol amd + BS')
%     tic
%     [R, ~, S] = chol(A);
%     tSetup = toc;
% 
%     tic
%     %C = S'*B;
%     B_cell = mat2tiles(S'*B, nv, BLOCK_SIZE);
%     parfor k = 1 : numel(B_cell)
%         x{k} = S * (R \ (R' \ B_cell{k}));
%     end
%     X = [x{:}];
%     %X = S * Z;
%     tSolve = toc;
%     tTotal = tSetup + tSolve;   
%     err = norm(A*X-B, 'fro') / norm(B, 'fro');
%     
%     TSetup(i, counter) = tSetup;
%     TSolve(i, counter) = tSolve;
%     TTotal(i, counter) = tTotal;
%     ERR(i, counter) = err;
%     counter = counter + 1;
%     if i == 1
%         names{end+1} = 'chol+amd+BS+parfor';
%     end
%     clear x B_cell X;
    
    %% --------------------------------------------------------------------
    % decomposition
    % ---------------------------------------------------------------------
%     disp('decomposition')
%     tic
%     dA = decomposition(A);
%     tDecompSetup = toc;
% 
%     tic
%     %%C = S'*B;
%     X = dA \ B;
%     %X = S * Z;
%     tDecompSolve = toc;
%     tCholTotal = tDecompSetup + tDecompSolve;   
%     cholErr = norm(A*X-B, 'fro') / norm(B, 'fro');
%     
%     TSetup(i, counter) = tDecompSetup;
%     TSolve(i, counter) = tDecompSolve;
%     TTotal(i, counter) = tCholTotal;
%     ERR(i, counter) = cholErr;
%     counter = counter + 1;
%     if i == 1
%         names{end+1} = 'decomposition';
%     end

    % ---------------------------------------------------------------------
    % ichol amd + BS
    % ---------------------------------------------------------------------
    disp('ichol amd + BS')
    tic
    Ri = ichol(Ap,struct('type','ict','droptol',1e-4,'michol','off'));
    tSetup = toc + tColamd;

    tic
    currentParms = spparms;
    spparms('autoamd',0);
    spparms('autommd',0);
    C = B(p, :);
    Z = Ri' \ (Ri \ C);
    X = Z(pt, :);
    spparms(currentParms);
    tSolve = toc;
    tTotal = tSetup + tSolve;
    err = norm(A*X-B, 'fro') / norm(B, 'fro');
    
    TSetup(i, counter) = tSetup;
    TSolve(i, counter) = tSolve;
    TTotal(i, counter) = tTotal;
    ERR(i, counter) = err;
    counter = counter + 1;
    if i == 1
        names{end+1} = 'ichol+amd+BS';
    end
    
    % ---------------------------------------------------------------------
    % ichol amd + parfor
    % ---------------------------------------------------------------------
    disp('ichol amd + parfor')
    tic
    Ri = ichol(Ap,struct('type','ict','droptol',1e-4,'michol','off'));
    tSetup = toc + tColamd;

    clear z;
    tic
    currentParms = spparms;
    spparms('autoamd',0);
    spparms('autommd',0);
    C = B(p, :);
    C_cell = mat2tiles(C, nv, BLOCK_SIZE);
    parfor k = 1:numel(C_cell)
        z{k} = Ri' \ (Ri \ C_cell{k});
    end
    X = [z{:}];
    X = X(pt, :);
    %Z = Ri' \ (Ri \ C);
    %X = Z(pt, :);
    spparms(currentParms);
    tSolve = toc;
    tTotal = tSetup + tSolve;
    err = norm(A*X-B, 'fro') / norm(B, 'fro');
        
    TSetup(i, counter) = tSetup;
    TSolve(i, counter) = tSolve;
    TTotal(i, counter) = tICholParforTotal;
    ERR(i, counter) = err;
    counter = counter + 1;
    if i == 1
        names{end+1} = 'ichol+amd+parfor+BS';
    end
    clear z C_cell;
    
    % ---------------------------------------------------------------------
    % ichol amd + pcg (tol=1e-6)
    % ---------------------------------------------------------------------
    disp('ichol amd + pcg (tol=1e-6)')
    tic
    Ri = ichol(Ap,struct('type','ict','droptol',DROPTOL,'michol','off'));
    tSetup = toc + tColamd;

    tic
    currentParms = spparms;
    spparms('autoamd',0);
    spparms('autommd',0);
    X = zeros(nv, K);
    for j = 1:K
        c = B(p, j);
        %z = Ri' \ (Ri \ c);
        [z, flag] = pcg(Ap, c, 1e-6, 100, Ri, Ri');
        X(:, j) = z;
    end
    X = X(pt, :);
    spparms(currentParms);
    tSolve = toc;
    tTotal = tSetup + tSolve;
    err = norm(A*X-B, 'fro') / norm(B, 'fro');
    
    TSetup(i, counter) = tSetup;
    TSolve(i, counter) = tSolve;
    TTotal(i, counter) = tTotal;
    ERR(i, counter) = err;
    counter = counter + 1;
    if i == 1
        names{end+1} = 'ichol+amd+pcg_e-6';
    end
    
    disp('ichol amd + pcg (tol=1e-4)')
    tic
    Ri = ichol(Ap,struct('type','ict','droptol',DROPTOL,'michol','off'));
    tSetup = toc + tColamd;

    tic
    currentParms = spparms;
    spparms('autoamd',0);
    spparms('autommd',0);
    X = zeros(nv, K);
    for j = 1:K
        c = B(p, j);
        %z = Ri' \ (Ri \ c);
        [z, flag] = pcg(Ap, c, 1e-4, 100, Ri, Ri');
        X(:, j) = z;
    end
    X = X(pt, :);
    spparms(currentParms);
    tSolve = toc;
    tTotal = tSetup + tSolve;
    err = norm(A*X-B, 'fro') / norm(B, 'fro');
    
    TSetup(i, counter) = tSetup;
    TSolve(i, counter) = tSolve;
    TTotal(i, counter) = tTotal;
    ERR(i, counter) = err;
    counter = counter + 1;
    if i == 1
        names{end+1} = 'ichol+amd+pcg_e-4';
    end
    
    disp('ichol amd + pcg (tol=1e-2)')
    tic
    Ri = ichol(Ap,struct('type','ict','droptol',DROPTOL,'michol','off'));
    tSetup = toc + tColamd;

    tic
    currentParms = spparms;
    spparms('autoamd',0);
    spparms('autommd',0);
    X = zeros(nv, K);
    for j = 1:K
        c = B(p, j);
        %z = Ri' \ (Ri \ c);
        [z, flag] = pcg(Ap, c, 1e-2, 100, Ri, Ri');
        X(:, j) = z;
    end
    X = X(pt, :);
    spparms(currentParms);
    tSolve = toc;
    tTotal = tSetup + tSolve;
    err = norm(A*X-B, 'fro') / norm(B, 'fro');
    
    TSetup(i, counter) = tSetup;
    TSolve(i, counter) = tSolve;
    TTotal(i, counter) = tTotal;
    ERR(i, counter) = err;
    counter = counter + 1;
    if i == 1
        names{end+1} = 'ichol+amd+pcg_e-2';
    end
    
    % ---------------------------------------------------------------------
    % ichol amd michol + BS
    % ---------------------------------------------------------------------
    disp('ichol amd + BS + michol')
    tic
    alpha = 1.0001;
    Ri = ichol(Ap,...
        struct('type','ict', ...
               'droptol',DROPTOL, ...
               'michol','on', ...
               'diagcomp', alpha));
    tSetup = toc + tColamd;

    tic
    currentParms = spparms;
    spparms('autoamd',0);
    spparms('autommd',0);
    C = B(p, :);
    Z = Ri' \ (Ri \ C);
    X = Z(pt, :);
    %X = Ri' \ (Ri \ B);
    spparms(currentParms);
    tSolve = toc;
    tTotal = tSetup + tSolve;
    err = norm(A*X-B, 'fro') / norm(B, 'fro');
    
    TSetup(i, counter) = tSetup;
    TSolve(i, counter) = tSolve;
    TTotal(i, counter) = tTotal;
    ERR(i, counter) = err;
    counter = counter + 1;
    if i == 1
        names{end+1} = 'ichol+michol+amd+BS';
    end
    
    % ---------------------------------------------------------------------
    % ichol amd michol + pcg
    % ---------------------------------------------------------------------
    disp('ichol amd + pcg + michol')
    tic
    alpha = 1.0001;
    Ri = ichol(Ap,...
        struct('type','ict',...
               'droptol',DROPTOL, ...
               'michol','on', ...
               'diagcomp', alpha));
    tSetup = toc + tColamd;

    tic
    currentParms = spparms;
    spparms('autoamd',0);
    spparms('autommd',0);
    X = zeros(nv, K);
    for j = 1:K
        c = B(p, j);
        %z = Ri' \ (Ri \ c);
        [z, flag] = pcg(Ap, c, 1e-8, 100, Ri, Ri');
        X(:, j) = z;
    end
    X = X(pt, :);
    %C = B(p, :);
    %Z = Ri' \ (Ri \ C);
    %X = Z(pt, :);
    spparms(currentParms);
    tSolve = toc;
    tTotal = tSetup + tSolve;
    err = norm(A*X-B, 'fro') / norm(B, 'fro');
    
    TSetup(i, counter) = tSetup;
    TSolve(i, counter) = tSolve;
    TTotal(i, counter) = tTotal;
    ERR(i, counter) = err;
    counter = counter + 1;
    if i == 1
        names{end+1} = 'ichol+michol+amd+pcg';
    end
    
    N_FACES(end+1) = m.nF;
end

%%
pltSet = {'chol+amd+BS', false; ...
    %'chol+amd+BS+parfor', false; ...
    'ichol+amd+BS', true; ...
    'ichol+amd+parfor+BS', false; ...
    'ichol+amd+pcg_e-6', true; ...    
    'ichol+amd+pcg_e-4', true; ...    
    'ichol+amd+pcg_e-2', true; ...    
    'ichol+michol+amd+BS', false; ...
    'ichol+michol+amd+pcg', false};
toPlot = containers.Map(pltSet(:,1), pltSet(:,2));
inds = cellfun(@(x) x, pltSet(:, 2));

plotting_defaults(30,30,2);
figure(1); 
subplot(221); hold on
[xx, I] = sort(N_FACES);
for i = 1:size(TSetup, 2)
    if ~toPlot(names{i})
        continue
    end
    yy = TSetup(I, i);
    plot(xx, yy, '--x')
end
hold off
title('Setup time')
legend(names{inds}, 'Location', 'northwest', 'Interpreter', 'none')

subplot(222); hold on
for i = 1:size(TSolve, 2)
    if ~toPlot(names{i})
        continue
    end
    yy = TSolve(I, i);
    plot(xx, yy, '--x');
end
hold off
title(sprintf('Solve time (%d rhs''s)', K))
legend(names{inds}, 'Location', 'northwest', 'Interpreter', 'none')

subplot(223); hold on
for i = 1:size(TTotal, 2)
    if ~toPlot(names{i})
        continue
    end
    yy = TTotal(I, i);
    plot(xx, yy, '--x');
end
hold off
title('Total time')
legend(names{inds}, 'Location', 'northwest', 'Interpreter', 'none')

subplot(224); hold on
for i = 1:size(ERR, 2)
    if ~toPlot(names{i})
        continue
    end
    yy = ERR(I, i);
    plot(xx, yy, '--x');
end
hold off
title('|AX - B| / |B|')
legend(names{inds}, 'Location', 'northwest', 'Interpreter', 'none')

if SAVE
    filename = fullfile(out_folder, 'time_solvers');
    print(filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
    save(filename, 'ERR', 'N_FACES', 'TTotal', 'TSetup', 'TSolve');
end


