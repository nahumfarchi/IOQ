%% Time inverting the laplacian using sparse, full, or gpu matrices.

fp = '../data/spot_s4.off';
disp(fp)

m = Mesh();
m.loadTM(fp);
nV = m.nV;
gd = gpuArray();

Lsp = -cotmatrix(m.V, m.F); % sparse
Lfl = full(Lsp); % full
Lsp_gpu = Lsp;
Lfl_gpu = single(Lfl);

%%
disp('------- CPU --------')

try
    disp(['spqr(Lsp: ', num2str(...
        timeit(@() spqr(Lsp)))])
catch ME
    disp('spqr(Lsp: FAILED')
    disp(getReport(ME))
end

try
    disp(['spinv(Lsp): ', num2str(...
        timeit(@() spinv(Lsp)))])
catch ME
    disp(['spinv(Lsp): FAILED'])
    disp(getReport(ME))
end
        
try
    disp(['inv(Lsp+1/nV)-1/nV: ', num2str(...
        timeit(@() inv(Lsp+1/nV)-1/nV))])
catch ME
    disp(['inv(Lsp+1/nV)-1/nV: FAILED'])
    disp(getReport(ME))
end

try
    disp(['inv(Lfl+1/nV)-1/nV: ', num2str(...
        timeit(@() inv(Lfl+1/nV)-1/nV))])
catch ME
    disp(['inv(Lfl+1/nV)-1/nV: FAILED'])
    disp(getReport(ME))
end

try
    disp(['pinv(Lfl): ', num2str(...
        timeit(@() pinv(Lfl)))])
catch ME
    disp(['pinv(Lfl): FAILED'])
    disp(getReport(ME))
end

try
    disp(['pinv(Lsp): ', num2str(...
        timeit(@() pinv(Lsp)))])
catch ME
    disp(['pinv(Lsp): FAILED'])
    disp(getReport(ME))
end

%%
disp('------- GPU --------')

try
    Lsp_gpu = gpuArray(Lsp_gpu);
    disp(['inv(Lsp_gpu+1/nV)-1/nV: ', num2str(...
        gputimeit(@() inv(Lsp_gpu+1/nV)-1/nV))])
catch ME
    disp(['inv(Lsp_gpu+1/nV)-1/nV: FAILED'])
    disp(getReport(ME))
end
Lsp_gpu = gather(Lsp_gpu);

try
    Lfl_gpu = gpuArray(Lfl_gpu);
    disp(['inv(Lfl_gpu+1/nV)-1/nV): ', num2str(...
        gputimeit(@() inv(Lfl_gpu+1/nV)-1/nV))])
catch ME
    disp(['inv(Lfl_gpu+1/nV)-1/nV): FAILED'])
    disp(getReport(ME))
end
Lfl_gpu = gather(Lfl_gpu);

try
    Lfl_gpu = gpuArray(Lfl_gpu);
    disp(['pinv(Lfl_gpu): ', num2str(...
        gputimeit(@() pinv(Lfl_gpu)))])
catch ME
    disp(['pinv(Lfl_gpu): FAILED'])
    disp(getReport(ME))
end
Lfl_gpu = gather(Lfl_gpu);

try
    Lsp_gpu = gpuArray(Lsp_gpu);
    disp(['Block inv (sparse), block size=8000: ', num2str(...
        gputimeit(@() block_inv_gpu(Lsp+1/nV, 8000)-1/nV))])
catch ME
    disp(['Block inv (sparse), block size=8000: FAILED'])
    disp(getReport(ME))
end
Lsp_gpu = gather(Lsp_gpu);

try
    Lfl_gpu = gpuArray(Lfl_gpu);
    disp(['Block inv (full), block size=8000: ', num2str(...
        gputimeit(@() ...
        block_inv_gpu(single(Lfl)+1/nV, 8000)-1/nV))])
catch ME
    disp(['Block inv (full), block size=8000: FAILED'])
    disp(getReport(ME))
end
Lfl_gpu = gather(Lfl_gpu);

% disp('------- CPU --------')
%     T = timeit(@() inv(Lsp+1/nV)-1/nV);
%     T1(ii) = T;
%     disp(['Sparse', num2str(T)])
%     
%     T = timeit(@() inv(Lfl+1/nV)-1/nV);
%     T2(ii) = T;
%     disp(['Full', num2str(T)])
% 
%     T = timeit(@() pinv(Lfl));
%     T3(ii) = T;
%     disp(['pinv (full)', num2str(T)])
% 
%     disp('------- GPU --------')
%     T = gputimeit(@() inv(Lsp_gpu+1/nV)-1/nV);
%     T4(ii) = T;
%     disp(['inv (sparse)', num2str(T)])
% 
%     T = gputimeit(@() inv(Lfl_gpu+1/nV)-1/nV);
%     T5(ii) = T;
%     disp(['inv (full)', num2str(T)])
% 
%     T = gputimeit(@() pinv(Lfl_gpu));
%     T6(ii) = T;
%     disp(['pinv (full)', num2str(T)])
% 
%     T = gputimeit(@() block_inv_gpu(Lsp_gpu+1/nV, 8000)-1/nV);
%     T7(ii) = T;
%     disp(['Block inv (sparse), block size=8000', num2str(T)])
% 
%     T = gputimeit(@() block_inv_gpu(Lfl_gpu+1/nV, 8000)-1/nV);
%     T8(ii) = T;
%     disp(['Block inv (full), block size=8000', num2str(T)])
