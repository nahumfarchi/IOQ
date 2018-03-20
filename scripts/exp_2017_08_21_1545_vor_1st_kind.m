function [] = exp_2017_08_21_1545_vor_1st_kind()
ME = [];
try
%MESHES = {'sphere_s0.off'};
ABOUT = 'VF using closest point in a lattice of Voronoi''s first kind\r\n\r\n';
OUT_FOLDER_NAME = 'voronoi';

%MESHES = {'sphere_s0.off', ...
%    'bumpy.off', ...
%    'round_cuber.off', ...
%    'rounded_cube_keenan.off', ...
%    'cow.off', ...
%    'bunny.off', ...
%    'bunny2.off', ...
%    'phands.off', ...
%    'torus_fat_r2.off'};

%MESHES = {'bunny2.off', 'bunny.off', 'rounded_cube_keenan.off', 'round_cuber.off'};
%MESHES = {'rounded_cube_keenan.off'};
%MESHES = {'round_cuber.off'};
%MESHES = {'sphere_s0.off'};
%MESHES = {'sphere_s0.off'};
MESHES = {'sphere_s0.off'};

%VIEW_ANGLE = [];

VERBOSE = true;
EPS = 1e-9;

theta0 = 0;

RUN_LIBIGL_MIQ = true;
PLOT = false;
SAVE = false;

if PLOT && SAVE
    OUT_FOLDER = create_time_stamped_folder(fullfile('..', 'results'), ...
        OUT_FOLDER_NAME, ...
        true, ...
        false);
    LOG = fopen(fullfile(OUT_FOLDER, 'log.txt'), 'w');
else
    LOG = -1;
end

log_and_print(LOG, '%s\r\n', mfilename('fullpath'));
log_and_print(LOG, ABOUT);
log_and_print(LOG, 'logid : %d\r\n', LOG);

k = 1;
for fname = MESHES
    log_and_print(LOG, 'Loading %s...\r\n', fname{:});
    p = find_data_folder();
    fp = fullfile(p, fname{:});

    mesh = Mesh();
    mesh.loadTM(fp);
    
    [d0, d1] = get_exterior_derivatives(mesh);
    Ad = get_gaussian_curvature(mesh);
    degree = 4;
    f0 = [1]; % Starting face
    v0 = mesh.V(mesh.F(f0, 2), :) - mesh.V(mesh.F(f0, 1), :);
    xi = 2 - 2*mesh.genus;

    log_and_print(LOG, 'degree : %d\r\n', degree);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % libigl MIQ (greedy rounding) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if RUN_LIBIGL_MIQ
        log_and_print(LOG, '\r\nSolving with libigl MIQ...\r\n');
        MIQ = nrosy_wrapper(fp, f0, v0, degree);
        log_and_print(LOG, '%s\r\n', MIQ.result);
        log_and_print(LOG, 'Status: %d\r\n', MIQ.status);
        log_and_print(LOG, 'E: %g\r\n', MIQ.E);
    end

    % |Mk-t|^2
    nV = mesh.nV;
    W = d0'*d0;
    [V,D] = eig(full(W));
    [D, sort_eigen] = sort(diag(D));
    V = V(:, sort_eigen);
    D = diag(D);

    D = D(1:nV, 1:nV);
    V = V(:, 1:nV);

    Dinv = diag(D);
    Dinv(Dinv < EPS) = inf;
    Dinv = 1 ./ Dinv;
    Dinv = diag(Dinv);
    Winv = sqrt(Dinv)*V';

    % |M*x-y|
    % |M*(x-z)|
    %M = V*sqrt(Dinv)*V';
    %b = (2/pi)*Ad;
    %t = (2/pi)*M*Ad;
    %m = 4*xi;
    
    % Resistence distance matrix
%     Wii = repmat(diag(Winv), 1, nV);
%     %G = -(Wii + Wii' - Winv - Winv');
%     %
%     G = (Wii + Wii' - 2*Winv);
%     G(1:nV+1:nV^2) = -sum(G, 2);
%     B = chol(G);
%     B = B(1:end-1,:);
%     %M = sqrt(G);

    B = Winv;
    
    z = (2/pi)*Ad;
    y = (2/pi)*B*Ad;
    c = 4*xi;
    
    %[X, b] = add_constraint(B, y, c);
    %k = closest_lpoint_vor(X, b);
    %k = closest_lpoint_vor(B, y);
    %u0 = MIQ.k;
    %u0(u0 > EPS) = 1;
    %u0(u0 < -EPS) = -1;
    %Z = B'*X;
    %D = inv(B');
    %tmp = B'*D;
    
    % a'*x=c
    a = (ones(1,nV-1) / (B(:, 1:end-1)))';
    tmp = norm(a);
    a = a / tmp;
    c = c / tmp;
    % Project y onto a'*x=c
    y_proj = y - a*(a'*(y - c*a));
    % Project B unto a'*x=c
    B_proj = B-repmat(a, 1, nV).*repmat(a'*B, nV, 1);
    %tmp = tmp(:, 1:rank(tmp));
    [x, uk] = closest_lpoint_vor(B_proj, y_proj);
    xx = x + c*a;
    u = tmp(:,1:end-2) \ xx;
    [~, I] = sort(u);
    k = zeros(nV, 1);
    k(I(end:-1:end-7)) = 0.25;
    
    %delta_ = 0.75;
    %tic
    %Ar = LLL_reduction(M(:,1:end-1), delta);
    %toc
    %
    %tic
    %k = babai(triu(Ar),t);
    %toc
    
    %k = LLL([M', ones(size(M',1),1)], [t', 4*xi]);
    %k = lattice_solve(M', t');
    %k = k(:);
    
%     % LLL_FP treats *rows* as basis vectors, while we treat *columns*.
%     M = LLL_FP(M')';
     %M = LLL_reduction(M);
%      M_reduced = CLLL(M);
%      M_B = M_reduced(:, 1);
%      M_I = M_reduced(:, 2:end);
%      c = 4*xi;
%      A = M_I - M_B*ones(1, numel(M_B)-1);
%      y = M_B*c - (M_reduced \ t);
%      [Q, R] = qr(A);
%      [Q2, R2] = qr(Q'*A); 
%      yy = R2\y;
%      k = babai(R2, yy);
%     k2 = babai2(A, y);
%     %[ksorted, inds] = sort(k);
%     %kmax = k(inds(end:-1:end-4));
%     %kmin = k(inds(1:4));
    
%     M_B = M(:, 1);
%     M_I = M(:, 2:end);
%     c = 4*xi;
%     A = M_I - M_B*ones(1, numel(M_B)-1);
%     y = M_B*c - t;
%     k = lattice_solve(A', y')';
%     k = [k; c-ones(1, numel(k))*k];
    
    %inds = find(abs(k) > EPS);
    %inds = inds(end:-1:end-7);
    %inds = inds(1:8);
    inds = find(abs(k) > EPS);
    S = [inds, k(inds)];
    LATTICE = TCODS(mesh, S, f0, theta0, degree, VERBOSE);
    fprintf('E = %g\n', LATTICE.E);
    
    LATTICE.title = {'Lattice lll', sprintf('$E_{MIQ} = %g$', LATTICE.E)};
    MIQ.title = {'MIQ ', sprintf('$E_{MIQ} = %g$', MIQ.E)};
    
    figure(1); clf(1)
    subplot(221)
    MeshVis.plot(mesh, ...
        'nrosy', LATTICE, ...
        'Title', LATTICE.title);
    subplot(222)
    MeshVis.plot(mesh, ...
        'nrosy', MIQ, ...
        'Title', MIQ.title);
end

% A hack since matlab does not have a finally clause
catch ME
end
% Close open resources
if LOG > 0
    fclose(LOG);
end
if ~isempty(ME)
    rethrow(ME);
end

    
function [X, b] = add_constraint(M, y, c)
    % x_I^{*} = argmin |X x - b|
    %           x \in Z^{n-1}
    % x^{*} = [x_I^{*}; ones(1, n-1) x_I^{*}]
    n = size(M, 2);
    M_B = M(:, 1);
    M_I = M(:, 2:end);
    X = M_I - M_B*ones(1,n-1);
    b = y - M_B*c;