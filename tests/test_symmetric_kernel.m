%% ========================================================================
%
%  ========================================================================
nv = 5;
L = magic(nv);
L = tril(L, -1);
L = gpuArray(single(L + L'));

Lp = inv(L+1/nv);
dl = diag(Lp);
R = bsxfun(@plus, dl, dl') - 2*Lp;
mask = tril(true(nv, nv), -1);
Rtril = R(mask);
b = single(zeros(nv, 1, 'gpuArray'));

opt.KernelEntry = 'reduce_cols_R_symmetric';
kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', opt.KernelEntry);
num_elem = nv^2;
kernel.ThreadBlockSize = [1, kernel.MaxThreadsPerBlock, 1];
kernel.GridSize = [1, ceil(nv / kernel.MaxThreadsPerBlock)];
out_min = single(zeros(1, nv, 'gpuArray'));
out_r = uint32(zeros(1, nv, 'gpuArray'));

[min_val1, i1] = min(R);
[out_min, out_r] = feval(kernel, out_min, out_r, Rtril, b, nv, nv, num_elem);
out_r = out_r + 1;
assert(norm(double(out_r) - i1) < 1e-10)
assert(norm(min_val1 - out_min) < 1e-6)


%%
%fp = '../data/bunny.off';
fp = '../data/bunnies/bunny_57k_faces.off';
m = Mesh(fp);
ne = m.nE; 
nv = m.nV; nf = m.nF; V = m.V; F = m.F;
d0 = get_exterior_derivatives(m);
L = d0'*d0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% True resistance
% print_header('gpu inv...');
% try
%     tic; Lp = inv(gpuArray(single(L+1/nv))) - 1/nv; wait(gd); toc
% catch
%     try
%         tic; Lp = block_inv_gpu(L+1/nv, 20000, true) - 1/nv; wait(gd); toc
%         
%     catch
%         tic; Lp = invChol_mex(full(L+1/nv)) - 1/nv; toc;
%     end
% end
% dl = diag(Lp);
%gputimeit(@() inv(gpuArray(L+1/nv)) - 1/nv)
%gputimeit(@() inv(gpuArray(single(L+1/nv))) - 1/nv) % faster


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Approx resistance

%% setup kernel
%kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', 'reduce_cols_Ztilde');
%kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', 'reduce_cols_R');
%nv = 4;
kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', 'reduce_cols_R_symmetric');
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

%% 
print_header('symmetric kernel');

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
elapsed1 = toc;

L = d0'*d0;

tic
t = 0.3;
L = L + (1/t)*speye(size(L, 1));

% create truncated Laplacian by removing positive off-diagonals
Lt = L - spdiags(diag(L), 0, size(L, 1), size(L, 2));
[i j v] = find(Lt);
v(v > 0) = 0;
Lt = sparse(i, j, v);
% add back excess diagonal!
Lt = Lt - spdiags(sum(Lt, 2) - 1/t, 0, size(L, 1), size(L, 2)); 

% setup preconditioner on the truncated Laplacian 
hsc_fun = hsc_setup(Lt, L);

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
    %zi = Lchol' \ (Lchol \ yi');
    zi = pcg(L, yi'/t, 1e-6, 100, hsc_fun, []);
    Ztilde(ii, :) = zi;
end
elapsed2 = toc;

L = d0'*d0;

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
tol = 1e-6; maxit = 100;
for ii = 1:k
    yi = Y(ii, :);
    zi = pcg(L, yi', tol, maxit, Lchol, Lchol');
    Ztilde(ii, :) = zi;
    %zi = Lchol' \ (Lchol \ yi');
    %Ztilde(ii, :) = zi';
    
end
elapsed3 = toc;

fprintf('Elapsed1, elapsed2 = %g, %g, %g\n', elapsed1, elapsed2, elapsed3);



%% kernel

Ztilde = gpuArray(single(Ztilde));
btilde = Lchol' \ (Lchol \ (2*(x-x0)));
btilde = gpuArray(single(btilde));

%Ztilde = gpuArray(single(magic(4)));
%btilde = gpuArray(single(zeros(4, 1)));
%btilde = gpuArray(single([10;2;3;4]));

%Rtilde = pdist(gpuArray(single(Ztilde'))).^2;
Rtilde = pdist(gpuArray(single(Ztilde')), 'squaredeuclidean');
wait(gd); fprintf('setup time : %g\n', toc);

tic
%tmp = single(squareform(Rtilde));
[out_min, out_r] = feval(kernel, out_min, out_r, Rtilde, btilde, nv, nv, num_elem);
out_r = out_r + 1;
[min_val1, j1] = min(out_min);
i1 = out_r(j1);
i1 = gather(i1);
wait(gd); fprintf('iter time : %g\n', toc);

%% bsx
[min_val, i2] = min(bsxfun(@plus, btilde, -btilde') + squareform(Rtilde));
[min_val, j2] = min(min_val);
i2 = i2(j2);

Ztilde = gather(Ztilde);
Rtilde = gather(Rtilde);
btilde = gather(btilde);

assert(i1 == i2, 'ii != i2')
assert(j1 == j2, 'j1 != j2')





