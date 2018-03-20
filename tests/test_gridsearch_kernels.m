fp = '../data/bunny.off';
%fp = '../data/bunnies/bunny_57k_faces.off';
m = Mesh(fp);
ne = m.nE; nv = m.nV; nf = m.nF; V = m.V; F = m.F;
d0 = get_exterior_derivatives(m);
L = d0'*d0;
disp('gpu inv...')
try
    tic; Lp = inv(gpuArray(single(L+1/nv))) - 1/nv; wait(gd); toc
catch
    try
        tic; Lp = block_inv_gpu(L+1/nv, 20000, true) - 1/nv; wait(gd); toc
        
    catch
        tic; Lp = invChol_mex(full(L+1/nv)) - 1/nv; toc;
    end
end
dl = diag(Lp);
%gputimeit(@() inv(gpuArray(L+1/nv)) - 1/nv)
%gputimeit(@() inv(gpuArray(single(L+1/nv))) - 1/nv) % faster

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% True resistance

%% setup kernel
kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', 'reduce_cols');
num_elem = nv^2;
kernel.ThreadBlockSize = [1, kernel.MaxThreadsPerBlock, 1];
kernel.GridSize = [1, ceil(nv / kernel.MaxThreadsPerBlock)];
out_min = single(zeros(1, nv, 'gpuArray'));
out_r = uint32(zeros(1, nv, 'gpuArray'));
gd = gpuDevice();

%% setup singularities
alpha_G = get_gaussian_curvature(m);
genus = round(1 - sum(alpha_G)/(4*pi));
x0 = (2/pi)*alpha_G;
c = round(abs(sum((x0)))); % = 4 xi

alpha_P = zeros(nv, 1);
n_pos_sing = (c + round(sum(x0))) / 2;
n_neg_sing = (c - round(sum(x0))) / 2;
inds_pos = randperm(nv, n_pos_sing);
inds_neg = randperm(nv, n_neg_sing);
alpha_P(inds_pos) = 1;
alpha_P(inds_neg) = -1;
x = alpha_P;

%% find min with kernel
Lp = gpuArray(single(Lp));
b = gpuArray(single(2*Lp*(x - x0)));
dl = gpuArray(single(dl));

disp('True R, kernel...')
tic
[out_min, out_r] = feval(kernel, out_min, out_r, Lp, dl, b, nv, nv, num_elem);
out_r = out_r + 1;
[min_val, j1] = min(out_min);
i1 = out_r(j1);
i1 = gather(i1);
wait(gd); toc

%% find min with bsx
disp('True R, bsx...')
tic
[min_val2, i2] = min(bsxfun(@plus, dl+b, (dl-b)') - 2*Lp);
[min_val2, j2] = min(min_val2);
i2 = i2(j2);
wait(gd); toc

assert(i1 == i2, 'i''s are not equal')
assert(j1 == j2, 'j''s are not equal')
assert(norm(min_val - min_val2) < 1e-6, 'min_val''s are not equal')

b = gather(b);
Lp = gather(Lp);
dl = gather(dl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Approx resistance

%% setup kernel
%kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', 'reduce_cols_Ztilde');
kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', 'reduce_cols_R');
num_elem = nv^2;
kernel.ThreadBlockSize = [1, kernel.MaxThreadsPerBlock, 1];
kernel.GridSize = [1, ceil(nv / kernel.MaxThreadsPerBlock)];
out_min = single(zeros(1, nv, 'gpuArray'));
out_r = uint32(zeros(1, nv, 'gpuArray'));
gd = gpuDevice();

%%
rng(112);

% setup singularities
alpha_G = get_gaussian_curvature(m);
genus = round(1 - sum(alpha_G)/(4*pi));
x0 = (2/pi)*alpha_G;
c = round(abs(sum((x0)))); % = 4 xi

alpha_P = zeros(nv, 1);
n_pos_sing = (c + round(sum(x0))) / 2;
n_neg_sing = (c - round(sum(x0))) / 2;
inds_pos = randperm(nv, n_pos_sing);
inds_neg = randperm(nv, n_neg_sing);
alpha_P(inds_pos) = 1;
alpha_P(inds_neg) = -1;
x = alpha_P;

%% Calculate Ztilde
print_header('Ztilde');
tic
Lchol = ichol(L, struct('type','ict','droptol',1e-4,'michol','off'));
eps = 0.5;
JLFac = 24;
k = round(JLFac * log(ne) / eps^2);
%k = round(JLFac * log(ne / eps^2));
% Approximate Z = Q d0 Lp
Q = 2*(rand(k, ne) > 0.5) - 1;
Q = (1/sqrt(k)) * Q;
Y = Q*d0;
Ztilde = zeros(k, nv);
for ii = 1:k
    yi = Y(ii, :);
    zi = Lchol' \ (Lchol \ yi');
    Ztilde(ii, :) = zi';
end
toc

%% find min with kernel
b = gpuArray(single(2*Lp*(x - x0)));
Ztilde = gpuArray(single(Ztilde));
btilde = Lchol' \ (Lchol \ (2*(x-x0)));
btilde = gpuArray(single(btilde));

disp('Approx R, kernel...')
tic
[out_min, out_r] = feval(kernel, out_min, out_r, Ztilde, k, btilde, nv, nv, num_elem);
out_r = out_r + 1;
[min_val3, j3] = min(out_min);
i3 = out_r(j3);
i3 = gather(i3);
wait(gd); toc

%% Find min with bsx
disp('Approx R, bsx...')
Ztilde = gather(Ztilde);
tic
Rtilde4 = sq_distance2(Ztilde, Ztilde);
[min_val4, i4] = min(bsxfun(@plus, btilde, -btilde') + Rtilde4);
[min_val4, j4] = min(min_val4);
i4 = i4(j4);
wait(gd); toc

assert(i3 == i4, 'i''s are not equal')
assert(j3 == j4, 'j''s are not equal')
assert(norm(min_val3 - min_val4) < 1e-6, 'min_val''s are not equal')

b = gather(b);
Ztilde = gather(Ztilde);
btilde = gather(btilde);
Rtilde4 = gather(Rtilde4);

%%
print_header('symmetric kernel');
Ztilde = gpuArray(single(Ztilde));
%btilde = Lchol' \ (Lchol \ (2*(x-x0)));
%btilde = gpuArray(single(btilde));

tic
Rtilde5 = pdist(Ztilde').^2;
%tmp = triu(Rtilde4, 1);
%tmp = tmp(tmp~=0);
%Rtilde5 = tmp;
[out_min, out_r] = feval(kernel, out_min, out_r, Rtilde5, btilde, nv, nv, num_elem);
out_r = out_r + 1;
[min_val, j5] = min(out_min);
i5 = out_r(j5);
i5 = gather(i5);
wait(gd); toc

Ztilde = gather(Ztilde);
Rtilde5 = gather(Rtilde5);
btilde = gather(btilde);

norm(squareform(Rtilde5) - Rtilde4, 'fro')





