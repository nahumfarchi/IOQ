
%% Timings and norm checks
BLOCK_SIZE = 8000;
gd = gpuDevice();
%fp = '../data/genus1_small/bucky_r_1000v.off';
fp = '../data/torus_s0.off';
m = Mesh(fp); nv = m.nV;
[d0, ~] = get_exterior_derivatives(m);
L = d0'*d0;

% with bucky180_r: 
%   T1 - too long
%   T2 - 662
%   T3 - 141.6716
%   T4 - out of mem
%   T5 - too long
%   T6 - 183.5232

%tic
%Lp_cpu = inv(L+1/nv) - 1/nv;
%T1 = toc 

tic
Lp_cpu_chol = invChol_mex(full(L+1/nv)) - 1/nv;
T2 = toc

%check_norm('Lp_cpu', 'Lp_cpu_chol', 'Log', -1, 'Tol', 1e-8);

tic
Lp_gpu_block = (block_inv_gpu(full(L+1/nv), BLOCK_SIZE) - 1/nv);
wait(gd); T3 = toc

check_norm('Lp_cpu_chol', 'Lp_gpu_block', 'Log', -1);

tic
Lp_gpu = gather(inv( gpuArray(L+1/nv) ) - 1/nv);
wait(gd); T4 = toc

check_norm('Lp_gpu', 'Lp_cpu_chol', 'Log', -1);

%tic
%Lp_gpu2 = gather(pinv(gpuArray(full(L))));
%wait(gd); T5 = toc

tic
Lp_gpu_block_no_mult = block_inv_gpu(full(L+1/nv), BLOCK_SIZE, false) - 1/nv;
wait(gd); T6 = toc

check_norm('Lp_gpu_block', 'Lp_gpu_block_no_mult', 'Log', -1);
check_norm('Lp_gpu_chol', 'Lp_gpu_block_no_mult', 'Log', -1);

%% more accurate timings (using timeit and gputimeit). Takes longer.
% TODO

