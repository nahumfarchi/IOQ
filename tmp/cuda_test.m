ker = parallel.gpu.CUDAKernel('add.ptx', 'add.cu', 'add2');
N = 8;
k.ThreadBlockSize = N;
in1 = ones(N,1,'gpuArray');
in2 = ones(N,1,'gpuArray');
result = feval(ker,in1,in2)
in1
in2
