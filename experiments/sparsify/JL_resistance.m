%% Check how close the approximated resistance distance is to the true resistance distance
meshname = 'torus_s0';
m = Mesh(['../../../data/', meshname, '.off']);
nv = m.nV; ne = m.nE;
eps = 0.5;
JLfac = 4.0;
%k = round(JLfac*log(nv / eps^2));
%k = round(JLfac * log(nv));
%k = nv;
%k = round(JLfac*log(2*ne)/eps^2);

N_REPS = 30;
%xx = (0:floor(nv/10):nv);
%if xx(end) ~= nv
%    xx(end+1) = nv;
%end
xx = ne;
E = zeros(length(xx), N_REPS);
T_R = zeros(length(xx), N_REPS);
T_R_tilde = zeros(length(xx), N_REPS);

[d0, ~] = get_exterior_derivatives(m);
L = d0' * d0;

for i = 1:length(xx)
    k = xx(i);
    %fprintf('nv : %d\neps : %g\n, k : %g\n', nv, eps, k);
    fprintf('%d // %d\n', i, length(xx));
    fprintf('k : %d\n', k);
    for j = 1:N_REPS
        fprintf('%d ', j);
        fprintf('R...\n');
        tic
        Lp = invChol_mex(full(L + 1/nv)) - 1/nv;
        dl = diag(Lp);
        R = bsxfun(@plus, dl, dl') - 2*Lp;
        T_R(i, j) = toc;

        % Calculate Z = Q * W^0.5 * B * Lp,
        % where Q is a random +-1/sqrt(k) matrix,
        %       W is a |E| x |E| diagonal weight matrix
        %       and B is d0
        fprintf('Rtilde...\n');
        tic
        Q = 2*(rand(k, ne) > 0.5) - 1;
        Q = (1 / sqrt(k)) * Q;
        %Q = randn(k, ne) * 1 / sqrt(k);

        Z = Q*d0*Lp;
        Rtilde = pdist2(Z', Z').^2;
        T_R_tilde(i, j) = toc;

        E(i, j) = norm(R - Rtilde);
    end
    fprintf('\n');
end

%%
figure(1)
%cmap = cbrewer('div','PuOr',100);
cmap = cbrewer('div','RdYlGn',3);
colormap(cmap);
X = E;
yy = mean(X, 2);
ee = std(X, 0, 2);
errorbar(xx(2:end), yy(2:end), ee(2:end));
title(sprintf('The mean of norm(R - Rtilde) over %d iterations', N_REPS));
print([meshname, '_JL_R_err'], '-dpng', '-r300');
saveas(gcf, [meshname, '_JL_R_err'], 'fig');


%%
figure(2)
X1 = T_R;
X2 = T_R_tilde;
yy1 = mean(X1, 2);
ee1 = std(X1, 0, 2);
yy2 = mean(X2, 2);
ee2 = std(X2, 0, 2);

%errorbar(xx, yy1, ee1, xx, yy2, ee2);
%legend('R', 'R tilde')
title(sprintf('Time to calculate R and R tilde (mean of %d runs)', N_REPS));
print([meshname, '_JL_R_time'], '-dpng', '-r300');
saveas(gcf, [meshname, '_JL_R_time'], 'fig');

% R_approx = zeros(nv, nv);
% for i = 1:nv
%     for j = 1:nv
%         R_approx(i, j) = norm(Z(:, i) - Z(:, j))^2;
%     end
% end
% toc
% 
% tic; 
% R_approx2 = pdist2(Z', Z'); 
% toc
% 
% norm(R_approx - R_approx2.^2)
% 
% norm(R - R_approx)
% 
% norm(R - R_approx2.^2)
% 
% cmp = [R(:), R_approx(:), R_approx2(:)];

%% Fast gridsearch
meshname = 'sphere_s4';
fp = ['../../../data/', meshname, '.off'];
m = Mesh(fp);
nv = m.nV; ne = m.nE;


disp('tcods_gsystem')
[~, Kf, d0, ~, H] = tcods_gsystem(m.V, m.F);

tic
disp('Lap...')
L = d0' * d0;
disp('ichol...')
Lchol = ichol(L,struct('type','ict','droptol',1e-04,'michol','off'));

disp('Q...')
eps = 0.5;
%JLfac = 4.0;
JLfac = 24;
k = round(JLfac * log(nv) / eps^2);
%k = nv;
rng(112)
Q = 2*(rand(k, ne) > 0.5) - 1;
Q = (1 / sqrt(k)) * Q;

disp('Ztilde...')
Y = Q*d0;
Ztilde = zeros(k, nv);
for ii = 1:k
    yi = Y(ii, :);
    zi = Lchol' \ ( Lchol \ yi' );
    Ztilde(ii, :) = zi';
end
%Z = Q*d0*Lp;
disp('Rtilde...')
Rtilde = pdist2(Ztilde', Ztilde').^2;

% EV = m.EVAdj;
% II = zeros(ne, 1);
% JJ = zeros(ne, 1);
% VV = zeros(ne, 1);
% for eid = 1:ne
%     eid
%     v1 = EV(eid, 1);
%     v2 = EV(eid, 2);
%     nz = norm(Ztilde(:, i) - Ztilde(:, j))^2;
%     II(eid) = v1;
%     JJ(eid) = v2;
%     VV(eid) = norm(Ztilde(:, i) - Ztilde(:, j))^2;
% end
% Rtilde = sparse([II; JJ], [JJ; II], [VV; VV], nv, nv, 2*ne);

% VV = m.VVAdj;
% coordinate = Ztilde';
% [n, d] = size(coordinate);
% resi = sparse(VV * diag(1:n));
% resj = sparse(diag(1:n) * VV);
% res = sparse(n, n);
% res(f) = arrayfun(@(x) (coordinate(resi(f), x) - coordinate(resj(f), x)) .^ 2, 1:d);

% gridsearch setup
disp('gridsearch setup...')
alpha_G = Kf(1:nv); beta_G = Kf(nv+1:end);
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
%b = 2*Lp*(x - x0);
min_val = inf;
E = [];

% gridsearch
for iter = 1:10
    disp([iter, min_val])
    %xtilde = L \ (2*(x-x0));
    xtilde = Lchol' \ ( Lchol \ (2*(x-x0)) );
    
    %[min_val, i] = min(Rtilde);
    %[min_val, j] = min(min_val);
    %i = i(j);
    
    %[min_val, i] = min(bsxfun(@plus, xtilde, -xtilde') + Rtilde);
    %[min_val, j] = min(min_val);
    %i = i(j);
    
    for eid = 1:ne
        v1 = EV(eid, 1);
        v2 = EV(eid, 2);
        val = xtilde(v1) - xtilde(v2) + Rtilde(v1, v2);
        if val < min_val
            min_val = val;
            i = v1;
            j = v2;
        end
        val = xtilde(v2) - xtilde(v1) + Rtilde(v2, v1);
        if val < min_val
            min_val = val;
            i = v1;
            j = v2;
        end
    end
    
    if abs(min_val(1)) < 1e-10
        disp('Minimization complete because optimality tolerance was reached.')
        disp(['m : ', num2str(min_val)])
        break
    end
    
    if i == j
        disp('Minimization complete because i==j');
        break
    end
    
    x(i) = x(i) + 1;
    x(j) = x(j) - 1;
    %b = b + 2*Lp(:,i) - 2*Lp(:,j);
    
    %E(iter) = gather((x-x0)'*Lp*(x-x0));
end
toc

%E3 = E;
res3 = TCODS(m, 'k', x, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
res3.miq_energy


%% Attempt to gridsearch using the approximated resistance matrix
meshname = 'sphere_s2';
fp = ['../../../data/', meshname, '.off'];
m = Mesh(fp);
nv = m.nV; ne = m.nE;
[~, Kf, d0, ~, H] = tcods_gsystem(m.V, m.F);

L = d0' * d0;
%tic
%Lp = invChol_mex(full(L + 1/nv)) - 1/nv;
%toc

%Lp = ichol(L);
%Lp = pinv(full(Lp));
%Lp = Lp' * Lp;

% Emiq1 : 17.9302
% Emiq2 :18.2692
% Emiq : 23.6917
%tic
%[Lp, ~] = chol(L+1/nv);
%Lp = inv(full(gpuArray(Lp)));
%Lp = Lp' * Lp;
%toc

%Lp = ichol(L);
%Emiq1 : 17.9302
%Emiq2 :18.171
%Emiq : 23.6917
tic
Lchol = ichol(L,struct('type','ict','droptol',1e-04,'michol','off'));
%Lp = single(full(gpuArray(Lp)));
%Lp = block_inv_gpu(full(Lp), 20000);
Lp = inv(Lchol);
Lp = Lp' * Lp;
toc

%tic
%Lp = inv(full(gpuArray(L)+1/nv)) - 1/nv;
%toc

dl = diag(Lp);
R = bsxfun(@plus, dl, dl') - 2*Lp;

eps = 0.5;
%JLfac = 4.0;
JLfac = 20;
k = round(JLfac * log(ne) / eps^2);
%k = nv;

rng(112)
Q = 2*(rand(k, ne) > 0.5) - 1;
Q = (1 / sqrt(k)) * Q;
%Q = randn(k, ne) * 1 / sqrt(k);

Z = Q*d0*Lp;
Rtilde = pdist2(Z', Z').^2;
fprintf('|R - Rtilde| = %g\n', norm(R-Rtilde));

alpha_G = Kf(1:nv); beta_G = Kf(nv+1:end);
x0 = (2/pi)*alpha_G;
c = round(abs(sum((x0)))); % = 4 xi

alpha_P = zeros(nv, 1);
n_pos_sing = (c + round(sum(x0))) / 2;
n_neg_sing = (c - round(sum(x0))) / 2;
inds_pos = randperm(nv, n_pos_sing);
inds_neg = randperm(nv, n_neg_sing);
alpha_P(inds_pos) = 1;
alpha_P(inds_neg) = -1;

% Regular gridsearch

x = alpha_P;
b = 2*Lp*(x - x0);
min_val = inf;
E = [];
tic
for iter = 1:100
    disp([iter, min_val])
    %[min_val, i] = min(bsxfun(@plus, b, -b') + Rtilde);
    %[min_val, j] = min(min_val);
    %i = i(j);
    [min_val, i] = min(bsxfun(@plus, dl+b, (dl-b)') - 2*Lp);
    [min_val, j] = min(min_val);
    i = i(j);

    if abs(min_val(1)) < 1e-10
        disp('Minimization complete because optimality tolerance was reached.')
        disp(['m : ', num2str(min_val)])
        break
    end
    
    if i == j
        disp('Minimization complete because i==j');
        break
    end
    
    x(i) = x(i) + 1;
    x(j) = x(j) - 1;
    b = b + 2*Lp(:,i) - 2*Lp(:,j);
    
    E(iter) = gather((x-x0)'*Lp*(x-x0));
end
toc

E1 = E;
res1 = TCODS(m, 'k', x, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);

% Gridsearch with approximated resistance and true b
x = alpha_P;
b = 2*Lp*(x - x0);
min_val = inf;
E = [];
tic
for iter = 1:100
    disp([iter, min_val])
    [min_val, i] = min(bsxfun(@plus, b, -b') + Rtilde);
    [min_val, j] = min(min_val);
    i = i(j);
    %[min_val, i] = min(bsxfun(@plus, dl+b, (dl-b)') - 2*Lp);
    %[min_val, j] = min(min_val);
    %i = i(j);

    if abs(min_val(1)) < 1e-10
        disp('Minimization complete because optimality tolerance was reached.')
        disp(['m : ', num2str(min_val)])
        break
    end
    
    if i == j
        disp('Minimization complete because i==j');
        break
    end
    
    x(i) = x(i) + 1;
    x(j) = x(j) - 1;
    b = b + 2*Lp(:,i) - 2*Lp(:,j);
    
    E(iter) = gather((x-x0)'*Lp*(x-x0));
end
toc

E2 = E;
res2 = TCODS(m, 'k', x, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);

% Gridsearch with approximated resistance and approximated 2(x-x0)Lp
x = alpha_P;
b = 2*Lp*(x - x0);
min_val = inf;
E = [];
tic
for iter = 1:100
    disp([iter, min_val])
    %xtilde = L \ (2*(x-x0));
    xtilde = Lchol' \ ( Lchol \ (2*(x-x0)) );
    %[min_val, i] = min(Rtilde);
    %[min_val, j] = min(min_val);
    %i = i(j);
    [min_val, i] = min(bsxfun(@plus, xtilde, -xtilde') + Rtilde);
    [min_val, j] = min(min_val);
    i = i(j);

    if abs(min_val(1)) < 1e-10
        disp('Minimization complete because optimality tolerance was reached.')
        disp(['m : ', num2str(min_val)])
        break
    end
    
    if i == j
        disp('Minimization complete because i==j');
        break
    end
    
    x(i) = x(i) + 1;
    x(j) = x(j) - 1;
    %b = b + 2*Lp(:,i) - 2*Lp(:,j);
    
    E(iter) = gather((x-x0)'*Lp*(x-x0));
end
toc

E3 = E;
res3 = TCODS(m, 'k', x, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);

figure
plot(1:length(E1), E1, 1:length(E2), E2, 1:length(E3), E3, '--')
legend('True R, true b', 'Approx R, true b', 'Approx R, approx b')
title('Gridsearch energy using true resistance and approximated resistance')

FACE0=1;THETA0=0;DEGREE=4;GVEC=[1,0,0];

[theta, p, kmiq, R, local_frames, frame_diffs, elapsed, pFixed] = ...
        NRoSy_mex_bin(fp, FACE0, GVEC, DEGREE);
Emiq = E_MIQ(m, theta, frame_diffs, p, DEGREE);

fprintf('Emiq1 : %g\nEmiq2 :%g\nEmiq3 : %g\nEmiq : %g\n', res1.miq_energy, res2.miq_energy, res3.miq_energy, Emiq);

%% Fast gridsearch by approximating Lp with ichol
meshname = 'bunny';
fp = ['../../../data/', meshname, '.off'];
m = Mesh(fp);
nv = m.nV; ne = m.nE;


disp('tcods_gsystem')
[~, Kf, d0, ~, H] = tcods_gsystem(m.V, m.F);

% s3 0.024
% s4 0.052
disp('Lap...')
tic
L = d0' * d0;
toc

disp('ichol...')
% s3 12.417
% s4 209.215
tic
Lchol = ichol(L,struct('type','ict','droptol',1e-04,'michol','off'));
%Lchol = chol(L)';
Lpapprox = zeros(nv, nv);
ei = zeros(nv, 1);
for vid = 1:nv
    ei(vid) = 1;
    Lpapprox(:, vid) = Lchol' \ (Lchol \ ei);
    ei(vid) = 0;
end
toc

% s3 19.726
% s4 1425.52
tic
Lp = invChol_mex(full(L+1/nv)) - 1/nv;
toc
norm(Lpapprox - Lp, 'fro');


% disp('Q...')
% s3 4.419
% s4 158.088
tic
eps = 0.5;
%JLfac = 4.0;
JLfac = 24;
k = round(JLfac * log(nv) / eps^2);
%k = nv;
rng(112)
Q = 2*(rand(k, ne) > 0.5) - 1;
Q = (1 / sqrt(k)) * Q;

disp('Ztilde...')
Y = Q*d0;
Ztilde = zeros(k, nv);
for ii = 1:k
    yi = Y(ii, :);
    zi = Lchol' \ ( Lchol \ yi' );
    Ztilde(ii, :) = zi';
end
%Z = Q*d0*Lp;
disp('Rtilde...')
%Rtilde = pdist2(Ztilde', Ztilde').^2;
Rtilde = sq_distance(Ztilde, Ztilde);

EV = m.EVAdj;
II = zeros(ne, 1);
JJ = zeros(ne, 1);
VV = zeros(ne, 1);
fprintf('ne = %g\n', ne);
for eid = 1:ne
    if mod(eid, 10000)==0, fprintf('%d ', eid); end
    v1 = EV(eid, 1);
    v2 = EV(eid, 2);
    zn = 1/norm(Ztilde(:, v1) - Ztilde(:, v2))^2;
    II(eid) = v1;
    JJ(eid) = v2;
    VV(eid) = zn;
    
    %Rtilde(v1, v2) = zn;
    %Rtilde(v2, v1) = zn;
end
Rtilde2 = sparse([II; JJ], [JJ; II], [VV; VV], nv, nv, 2*ne);

toc

% EV = m.EVAdj;
% II = zeros(ne, 1);
% JJ = zeros(ne, 1);
% VV = zeros(ne, 1);
% for eid = 1:ne
%     eid
%     v1 = EV(eid, 1);
%     v2 = EV(eid, 2);
%     nz = norm(Ztilde(:, i) - Ztilde(:, j))^2;
%     II(eid) = v1;
%     JJ(eid) = v2;
%     VV(eid) = norm(Ztilde(:, i) - Ztilde(:, j))^2;
% end
% Rtilde = sparse([II; JJ], [JJ; II], [VV; VV], nv, nv, 2*ne);

% VV = m.VVAdj;
% coordinate = Ztilde';
% [n, d] = size(coordinate);
% resi = sparse(VV * diag(1:n));
% resj = sparse(diag(1:n) * VV);
% res = sparse(n, n);
% res(f) = arrayfun(@(x) (coordinate(resi(f), x) - coordinate(resj(f), x)) .^ 2, 1:d);

% Regular gridsearch

x = alpha_P;
b = 2*Lp*(x - x0);
dl = diag(Lp);
min_val = inf;
E = [];
tic
for iter = 1:100
    disp([iter, min_val])
    %[min_val, i] = min(bsxfun(@plus, b, -b') + Rtilde);
    %[min_val, j] = min(min_val);
    %i = i(j);
    [min_val, i] = min(bsxfun(@plus, dl+b, (dl-b)') - 2*Lp);
    [min_val, j] = min(min_val);
    i = i(j);

    if abs(min_val(1)) < 1e-10
        disp('Minimization complete because optimality tolerance was reached.')
        disp(['m : ', num2str(min_val)])
        break
    end
    
    if i == j
        disp('Minimization complete because i==j');
        break
    end
    
    x(i) = x(i) + 1;
    x(j) = x(j) - 1;
    b = b + 2*Lp(:,i) - 2*Lp(:,j);
    
    E(iter) = gather((x-x0)'*Lp*(x-x0));
end
toc

E1 = E;
res1 = TCODS(m, 'k', x, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);

% gridsearch setup
disp('gridsearch setup...')
alpha_G = Kf(1:nv); beta_G = Kf(nv+1:end);
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
%b = 2*Lp*(x - x0);
min_val = inf;
E = [];

% gridsearch
EV = m.EVAdj;
for iter = 1:30
    disp([iter, min_val])
    %xtilde = L \ (2*(x-x0));
    tic
    xtilde = Lchol' \ ( Lchol \ (2*(x-x0)) );
       
    %[min_val, i] = min(Rtilde);
    %[min_val, j] = min(min_val);
    %i = i(j);
    
    [min_val, i] = min(bsxfun(@plus, xtilde, -xtilde') + Rtilde);
    [min_val, j] = min(min_val);
    i = i(j);
    
%     tic
%     for eid = 1:ne
%         v1 = EV(eid, 1);
%         v2 = EV(eid, 2);
%         val = xtilde(v1) - xtilde(v2) + Rtilde(v1, v2);
%         if val < min_val
%             min_val = val;
%             i = v1;
%             j = v2;
%         end
%         val = xtilde(v2) - xtilde(v1) + Rtilde(v2, v1);
%         if val < min_val
%             min_val = val;
%             i = v1;
%             j = v2;
%         end
%     end
%     toc
    
    if abs(min_val(1)) < 1e-10
        disp('Minimization complete because optimality tolerance was reached.')
        disp(['m : ', num2str(min_val)])
        break
    end
    
    if i == j
        disp('Minimization complete because i==j');
        break
    end
    
    x(i) = x(i) + 1;
    x(j) = x(j) - 1;
    %b = b + 2*Lp(:,i) - 2*Lp(:,j);
    
    E(iter) = gather((x-x0)'*Lp*(x-x0));
end

E3 = E;
res3 = TCODS(m, 'k', x, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
res3.miq_energy

FACE0=1;THETA0=0;DEGREE=4;GVEC=[1,0,0];
[theta, p, kmiq, R, local_frames, frame_diffs, elapsed, pFixed] = ...
        NRoSy_mex_bin(fp, FACE0, GVEC, DEGREE);
Emiq = E_MIQ(m, theta, frame_diffs, p, DEGREE);

figure
plot(1:length(E1), E1*(pi/2)^2, 1:length(E3), E3*(pi/2)^2, 1:length(E3), ones(length(E3), 1) * Emiq, '--')
legend('true R', 'approx R', 'MIQ (final energy)')
title('Gridsearch energy (with respect to Lp) with true R and approx R')
print('gridsearch_E_R_vs_Rtilde', '-dpng', '-r300')


