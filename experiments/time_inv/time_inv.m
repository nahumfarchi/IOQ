files = {'../../../data/sphere_s0.off', ...
    '../../../data/sphere_s1.off', ...
    '../../../data/sphere_s2.off', ...
    '../../../data/sphere_s3.off'};
lgds = {'invChol_mex', 'chol'};
n_files = numel(files);
n_experiments = 2;
T = zeros(n_files, n_experiments);
ACC = zeros(n_files, n_experiments);
for i = 1:n_files
    j = 1;
    fp = files{i};
    m = Mesh(fp);
    nv = m.nV;
    [d0, ~] = get_exterior_derivatives(m);
    L = d0' * d0;
    
    disp('chol mex')
    tic
    Lp_chol_mex = invChol_mex(L+1/nv) - 1/nv;
    elapsed = toc
    T(i, j) = elapsed;
    %ACC(i, j) = norm(speye(nv) - Lp_chol_mex*L)
    clear Lp_chol_mex
    
%     disp('lu')
%     j = j + 1;
%     tic
%     [lower, upper] = lu(L+1/nv);
%     Lp_lu = inv(upper) * inv (lower);
%     elapsed = toc
%     T(i, j) = elapsed;
%     clear Lp_lu

    disp('chol')
    j = j + 1;
    tic
    R = chol(L+1/nv);
    Lp_chol = inv(R) * inv(R');
    elapsed = toc
    T(i, j) = elapsed;
    clear R
    clear Lp_chol
    
end

figure
hold on
for j = 1:n_experiments
    plot(T(:, j))
end
hold off
legend(lgds)