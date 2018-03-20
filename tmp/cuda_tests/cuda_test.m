k = parallel.gpu.CUDAKernel('test.ptx', 'test.cu', 'add_mat');
N = 32;
k.ThreadBlockSize = [N, N, 1];
Lp = ones(N, N, 'gpuArray');
%in2 = ones(N, N, 'gpuArray');
b = gpuArray(1:N)';
dlp = gpuArray((1:N)+4)';

tic
result1 = gather(feval(k, Lp, dlp, b, N, N, N^2));
t1 = toc()

tic
result2 = gather(bsxfun(@plus, dlp+b, (dlp-b)') - 2*Lp);
t2 = toc()

check_norm('result1', 'result2');