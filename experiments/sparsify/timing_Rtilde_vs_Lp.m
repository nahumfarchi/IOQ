%%
gd = gpuDevice();
fp = '../../../data/ashish_nob/buste.off';
m = Mesh(fp); nv = m.nV; ne = m.nE;

%%
disp('d0, L')
tic
d0 = get_exterior_derivatives(m);
L = d0' * d0;
toc

disp('Lp gpu');
tic
%inv(gpuArray(single(L+1/nv))) - 1/nv;
block_inv_gpu(single(full(L+1/nv)), 20000) - 1/nv;
wait(gd); toc

% disp('chol')
% tic
% chol(L);
% toc
% 
% disp('gpu chol')
% tic
% chol(gpuArray(single(full(L))));
% wait(gd); toc
% 
% disp('ichol')
% tic
% ichol(L,struct('type','ict','droptol',1e-04,'michol','off'));
% toc

disp('ichol')
tic
Lchol = ichol(L,struct('type','ict','droptol',1e-04,'michol','off'));
%Lchol = gpuArray(Lchol);
eps = 0.5;
%JLfac = 4.0;
JLfac = 24;
k = round(JLfac * log(nv) / eps^2);
%k = nv;
rng(112)
Q = 2*(rand(k, ne) > 0.5) - 1;
Q = (1 / sqrt(k)) * Q;

%disp('Ztilde...')
Y = Q*d0;
Ztilde = zeros(k, nv);
fprintf('k = %g\n', k);
for ii = 1:k
    if mod(ii,100)==0, fprintf('%d ', ii); end
    yi = Y(ii, :);
    zi = Lchol' \ ( Lchol \ yi' );
    Ztilde(ii, :) = zi';
end
fprintf('\n');

%Z = Q*d0*Lp;
%disp('Rtilde...')
%Rtilde = pdist2(Ztilde', Ztilde').^2;
%Ztilde = gpuArray(Ztilde);
Rtilde = sq_distance(Ztilde, Ztilde); % this seems to be the bottleneck

%Rtilde = sparse(nv, nv);
% EV = m.EVAdj;
% II = zeros(ne, 1);
% JJ = zeros(ne, 1);
% VV = zeros(ne, 1);
% fprintf('ne = %g\n', ne);
% for eid = 1:ne
%     if mod(eid, 10000)==0, fprintf('%d ', eid); end
%     v1 = EV(eid, 1);
%     v2 = EV(eid, 2);
%     zn = norm(Ztilde(:, v1) - Ztilde(:, v2))^2;
%     II(eid) = v1;
%     JJ(eid) = v2;
%     VV(eid) = zn;
%     
%     %Rtilde(v1, v2) = zn;
%     %Rtilde(v2, v1) = zn;
% end
% Rtilde = sparse([II; JJ], [JJ; II], [VV, VV], nv, nv, 2*ne);
% fprintf('\n')
toc
clear Lchol
clear Ztilde;