%%
%eps = 0.5;
%n = 256*100;
%m = floor(24*log(n)/eps^2);
%m = 256*10;

n = 300*256;
m = 4*256;
fprintf('n, m = %d, %d\n', n, m);

%A = magic(m);
A = rand(n, m);
A = A  ./ norm(A, 'fro');
%A = gpuArray(single(A));
A = gpuArray(A);
gd = gpuDevice();

% matlab
print_header('pdist');
tic
D1 = pdist(A).^2;
wait(gd); toc
D1 = gather(D1);
D1 = squareform(D1);
%fprintf('pdist time : %.4g\n', gputimeit(@() pdist2(A, A).^2));

% bsx fun
print_header('sq_distance2');
try
    tic
    D2 = sq_distance2(A', A');
    wait(gd); toc
    
    norm(D1 - D2, 'fro')
    D2 = gather(D2);
catch
    disp('sq_distance2 failed')
end

print_header('sq_distance');
try
    tic
    D22 = sq_distance(A');
    wait(gd); toc
    
    norm(D1 - D22, 'fro')
    D22 = gather(D22);
catch
    disp('sq_distance failed')
end
%fprintf('sq_distance time : %.4g\n', gputimeit(@() sq_distance(A', A')));

D1 = gather(D1);


%% gpuPdist1
%tic
%kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', 'gpuPdist1');
%kernel.ThreadBlockSize = 256;
%kernel.GridSize = ceil(n / 256);
%D3 = gpuArray(zeros(n, n, 'single'));
%D3 = feval(kernel, D3, A, n, m);
%wait(gd); toc

%% kernel (row major kernel)
A = A';
%[m, n] = deal(n, m);
SZ = 32;
tic
kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', 'matrix_euclidean_distance_kernel_fast');
%num_elem = nv^2;
kernel.ThreadBlockSize = [SZ, SZ];
kernel.GridSize = [ceil(n / SZ), ceil(n / SZ)];
%out_min = single(zeros(1, nv, 'gpuArray'));
%out_r = uint32(zeros(1, nv, 'gpuArray'));
%D4 = gpuArray(zeros(n, n, 'single'));
D4 = gpuArray(zeros(n, n, 'single'));
D4 = feval(kernel, D4, A, n, m);
wait(gd); toc

D4 = gather(D4);

norm(D1 - D4, 'fro')

%[m, n] = deal(n, m);
A = A';
%clear D4

%clear A D1 D2 D3;

%% kernel col reduce (col major kernel)
tic
nelem = n * m;
SZ = 16; % threads per block
bpg = ceil( nelem / SZ ); % blocks per grid
kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', 'gpuPdist_col');
kernel.ThreadBlockSize = [SZ, SZ];
kernel.GridSize = [bpg, bpg];
D5 = gpuArray(zeros(n, n, 'single'));
D5 = feval(kernel, D5, A, n, m);
wait(gd); toc

D5 = gather(D5);

norm(D1 - D5, 'fro')

%% gpuPdist2
tic
A = A';
sz = 32;
kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', 'gpuPdist2');
kernel.ThreadBlockSize = [sz, sz];
kernel.GridSize = [ceil(n / sz), ceil(n / sz)];
D6 = gpuArray(zeros(n, n, 'single'));
D6 = feval(kernel, D6, A, n, m);
A = A';
wait(gd); toc

D6 = gather(D6);
norm(D1 - D6, 'fro')