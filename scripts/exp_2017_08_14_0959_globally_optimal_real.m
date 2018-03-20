function [] = exp_2017_08_14_0959_globally_optimal_real()
ME = [];
try
%MESHES = {'sphere_s0.off'};
ABOUT = 'Globally optimal, real with 1 constraint\r\n\r\n';
OUT_FOLDER_NAME = 'GO_real_with_1_constraint_N1';

MESHES = {'sphere_s0.off', ...
    'bumpy.off', ...
    'round_cuber.off', ...
    'rounded_cube_keenan.off', ...
    'cow.off', ...
    'bunny.off', ...
    'bunny2.off', ...
    'phands.off', ...
    'torus_fat_r2.off'};

%MESHES = {'bunny2.off', 'bunny.off', 'rounded_cube_keenan.off', 'round_cuber.off'};
%MESHES = {'rounded_cube_keenan.off'};
%MESHES = {'round_cuber.off'};
%MESHES = {'sphere_s0.off'};
%MESHES = {'sphere_s0.off'};

VIEW_ANGLE = [-11, 16];

VERBOSE = true;
EPS = 1e-10;

theta0 = 0;

RUN_LIBIGL_MIQ = true;
PLOT = true;
SAVE = true;

if PLOT && SAVE
    OUT_FOLDER = create_time_stamped_folder(fullfile('..', 'results'), ...
        OUT_FOLDER_NAME, ...
        true, ...
        false);
    LOG = fopen(fullfile(OUT_FOLDER, 'log.txt'), 'w');
else
    OUT_FOLDER = [];
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
    degree = 1;
    f0 = [1]; % Starting face
    v0 = mesh.V(mesh.F(f0, 2), :) - mesh.V(mesh.F(f0, 1), :);
    v0 = normalize_rows(v0);

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
    
       
    %%%%%%%%%%%%%%%%%%%%
    % Globally optimal %
    %%%%%%%%%%%%%%%%%%%%
    log_and_print(LOG, '\r\nBuilding GO system...\r\n');
    EF = mesh.EFAdj;
    [local_frames, diffs] = create_local_frames(mesh);
    r = exp(1i*degree*diffs);
    count = 1;
    
    I = [];
    J = [];
    S = [];
    nF = mesh.nF;
    % Construct smooth energy system E_s = norm(A*x)^2
    % where
    %   u_f = x(f)+1i*x(nF+f)
    %   u_g = x(g)+1i*x(nF+g)
    for eid = 1:mesh.nE
        if mesh.isBoundaryEdge(eid)
            error('Boundary edges not implemented')
        end
        
        f = EF(eid, 1);
        g = EF(eid, 2);
        
        % See notes/2017-08-02 globally optimal.pdf
        I(end+1) = count;
        J(end+1) = f;
        S(end+1) = cos(degree*diffs(eid));
        
        I(end+1) = count;
        J(end+1) = nF+f;
        S(end+1) = -sin(degree*diffs(eid));
        
        I(end+1) = count;
        J(end+1) = g;
        S(end+1) = -1;
        
        I(end+1) = count+1;
        J(end+1) = f;
        S(end+1) = sin(degree*diffs(eid));
        
        I(end+1) = count+1;
        J(end+1) = nF+f;
        S(end+1) = cos(degree*diffs(eid));
        
        I(end+1) = count+1;
        J(end+1) = nF+g;
        S(end+1) = -1;
        
        count = count + 2;
    end
    
    % Constraint energy E_l = -<phi,psi_0>
    k = 1;
    for fid = f0
        v = v0(k, :);
        vx = dot(v, local_frames(fid, :));
        vy = dot(v, local_frames(nF+fid, :));
        
        I(end+1) = count;
        J(end+1) = fid;
        S(end+1) = vx;
        
        I(end+1) = count;
        J(end+1) = nF+fid;
        S(end+1) = vy;
        
        count = count + 1;
        k = k + 1;
    end
    
    assert(2*mesh.nE + 1 == count-1);
    n = 2*mesh.nE + 1;
    m = 2*mesh.nF;
    X = sparse(I, J, S, n, m);
    %b = zeros(n);
    
    
    %[V, D] = eig(full(A'*A));
    
    %lam = 1;
    %lam = 0.1;
    %lam = 0.0227;
    %lam = -0.1;
    %lam = -1;
    %lam = -10;
    %lam = -100;
    %lam = -1000;
    %[tmp1, tmp2] = eig(full(A'*A));
    %lambdas = linspace(-10, tmp2(1,1), 10);
    
    %for lam = lambdas
    
    
    %[V, D] = eig((A-lam*eye(size(A)))'*(A-lam*eye(size(A))));
    %x = V(:, 1);
    %x = (A(1:end-1,:)'*A(1:end-1,:)) \ A(end,:)';
    %x = (A'*A - lam*eye(size(A'*A))) \ A(end,:)';
    %x = (A(1:end-1,:) - lam*eye(size(A(1:end-1,:)))) \ A(end,:)';
    %x = A(1:end-1,:)'*A(1:end-1,:)-lam*eye(size(A(1:end-1,:)'*A(1:end-1,:))) \ A(end,:)';
    
    %[V,D] = eig(full(A'*A)-lam*eye(size(A'*A)));
    %x = V(:, 1);
    
    log_and_print(LOG, '\r\nSolving...\r\n');
    %X = A'*A;
    A_s = X(1:end-1,:)'*X(1:end-1,:);
    %A_s = X'*X;
    phi0 = X(end,:)';
    
    %[V,D] = eigs(A_s, 6, 'SM');
    %[D,inds] = sort(diag(D));
    %V = V(:, inds);
    %D = diag(D);
    
    %[V,D] = eig(full(A_s));
    [V,D] = eigs(A_s, 6, 'sm');
    
    %phi0(1,1) = V(1,1);
    %phi0(nF+1,1) = V(nF+1,1);
    %phi0 = phi0 / norm(phi0);
    %v0 = local_frames(1,:)*phi0(1,1)+local_frames(1+nF,:)*phi0(1+nF,1);
    
    %lam = D(1,1);
    sz = 1;
    lambdas = linspace(-10, D(1,1)-1e-10, sz);
    %lambdas = D(1,1) - EPS;
    k = 1;
    GO = {};
    GO_proj = {};
    DD = diag(D) - D(1,1);
    n_lam1 = length(find(abs(DD)<EPS));
    %GO_TC = {};
    for lamt = lambdas
        log_and_print(LOG, '\r\nlam = %g\r\n', lamt);
        %x = (A(1:end-1,:)'*A(1:end-1,:)-lam*eye(size(A(1:end-1,:)'*A(1:end-1,:)))) \ A(end,:)';
        %x = (X'*X - lam*eye(size(X'*X))) \ X(end,:)';
        x = (A_s - lamt*speye(size(A_s))) \ phi0;
        check_norm('(A_s-lamt*speye(size(A_s)))*x', 'phi0');
    
        x = x / norm(x);
        U = V(:, 1:n_lam1);
        z = U*(U'*x);
        %norm(x-x_proj)
        %x = x_proj;
        
        check_norm('A_s*z', 'lamt*z'); % (a)
        
        u = zeros(mesh.nF, 1);
        u_proj = zeros(mesh.nF, 1);
        for fid = 1:mesh.nF
            a = x(fid);
            b = x(nF+fid);
            a_proj = z(fid);
            b_proj = z(nF+fid);
            %x(fid) = x(fid) / sqrt(a^2+b^2);
            %x(nF+fid) = x(nF+fid) / sqrt(a^2+b^2);
            u(fid) = a+1i*b;
            u_proj(fid) = a_proj + 1i*b_proj;
        end
        
        GO{k}.u = u;
        GO{k}.t = 1 / (1+norm((A_s - lamt*speye(size(A_s)))*phi0));
        %GO{k}.E_GO = (1-t)*norm(A_s*x) - t*dot(phi0, x);
        %GO{k}.E_GO = E_GO(mesh, u./abs(u), degree);
        GO{k}.E_GO = norm(A_s*x)^2;
        GO{k}.theta = angle(GO{k}.u) / degree;
        GO{k}.ffield = angles_to_ffield(GO{k}.theta, local_frames, degree);
        %GO{k}.ffield = repmat(abs(u), degree, 3) .* GO{k}.ffield;
        GO{k}.degree = degree;
        GO{k}.p = zeros(mesh.nE, 1);
        for eid = 1:mesh.nE
            fi = EF(eid, 1);
            fj = EF(eid, 2);
            GO{k}.p(eid) = round(-degree/(2*pi) * (GO{k}.theta(fi) + diffs(eid) - GO{k}.theta(fj)));
        end
        GO{k}.E_MIQ = E_MIQ(mesh, GO{k}.theta, diffs, GO{k}.p, degree);
        [GO{k}.x, GO{k}.k] = MIQ_to_TCoDS(GO{k}.theta, GO{k}.p, diffs, MIQ.Kg, d0, d1, degree);
        inds = find(abs(GO{k}.k) > EPS);
        GO{k}.S = [inds, GO{k}.k(inds)/degree];
    
        GO_TC = TCODS(mesh, GO{k}.S, f0, theta0, degree);
        
        GO{k}.title = {'GO (backslash)', ...
            sprintf('$E_{MIQ} = %g$', GO{k}.E_MIQ), ...
            sprintf('$E_{GO} = %g$', GO{k}.E_GO), ...
            sprintf('$\\lambda_{t}=%g$', lamt)};
        %GO_TC.title = {'TC (GO)', sprintf('$E_{MIQ} = %g$', GO_TC.E), sprintf('$E_{GO} = %g$', GO_TC.E_GO)};
        
        GO_proj{k}.u = u_proj;
        GO_proj{k}.t = 1 / (1+norm((A_s - lamt*speye(size(A_s)))*phi0));
        %GO_proj{k}.E_GO = (1-t)*norm(A_s*x) - t*dot(phi0, x);
        %GO_proj{k}.E_GO = E_GO(mesh, u./abs(u), degree);
        GO_proj{k}.E_GO = norm(A_s*x)^2;
        GO_proj{k}.theta = angle(GO_proj{k}.u) / degree;
        GO_proj{k}.ffield = angles_to_ffield(GO_proj{k}.theta, local_frames, degree);
        %GO_proj{k}.ffield = repmat(abs(u), degree, 3) .* GO_proj{k}.ffield;
        GO_proj{k}.degree = degree;
        GO_proj{k}.p = zeros(mesh.nE, 1);
        for eid = 1:mesh.nE
            fi = EF(eid, 1);
            fj = EF(eid, 2);
            GO_proj{k}.p(eid) = round(-degree/(2*pi) * (GO_proj{k}.theta(fi) + diffs(eid) - GO_proj{k}.theta(fj)));
        end
        GO_proj{k}.E_MIQ = E_MIQ(mesh, GO_proj{k}.theta, diffs, GO_proj{k}.p, degree);
        [GO_proj{k}.x, GO_proj{k}.k] = MIQ_to_TCoDS(GO_proj{k}.theta, GO_proj{k}.p, diffs, MIQ.Kg, d0, d1, degree);
        inds = find(abs(GO_proj{k}.k) > EPS);
        GO_proj{k}.S = [inds, GO_proj{k}.k(inds)/degree];
    
        %GO_TC = TCODS(mesh, GO_proj{k}.S, f0, theta0, degree);
        
        GO_proj{k}.title = {'z', ...
            sprintf('$E_{MIQ} = %g$', GO_proj{k}.E_MIQ), ...
            sprintf('$E_{GO} = %g$', GO_proj{k}.E_GO), ...
            sprintf('$\\lambda_{t}=%g$', lamt)};
        %GO_TC.title = {'TC (GO_proj)', sprintf('$E_{MIQ} = %g$', GO_TC.E), sprintf('$E_{GO} = %g$', GO_TC.E_GO)};
        
        k = k + 1;
    end
       
    % Calculate GO energy
    MIQ.u = langles_to_complex(MIQ.theta);
    tmp = [real(MIQ.u); imag(MIQ.u)];
    MIQ.E_GO = tmp'*A_s*tmp;
    %MIQ.E_GO = norm(X*[real(MIQ.u); imag(MIQ.u)]-b)^2;   
    GO_TC.u = langles_to_complex(GO_TC.theta);
    tmp = [real(GO_TC.u); imag(GO_TC.u)];
    GO_TC.E_GO = tmp'*A_s*tmp;
    %GO_TC.E_GO = norm(A*[real(GO_TC.u); imag(GO_TC.u)]-b)^2;
    
    % Plot
    MIQ.title = {'MIQ', ...
        sprintf('$E_{MIQ} = %g$', MIQ.E), ...
        sprintf('$E_{GO} = %g$', MIQ.E_GO)};
    GO_TC.title = {'TC (GO)', ...
        sprintf('$E_{MIQ} = %g$', GO_TC.E), ...
        sprintf('$E_{GO} = %g$', GO_TC.E_GO)};
    
    if PLOT
        MeshVis.wfigs(fname{:}, ...
            mesh, ...
            'OutFolder', OUT_FOLDER, ...
            'Nrosy', {MIQ, GO{:}, GO_proj{:}, GO_TC}, ...
            'ConstrainedFaces', f0, ...
            'ConstraintVectors', v0, ...
            'View', VIEW_ANGLE, ...
            'Save', SAVE);
    end
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

    function E_ = E_GO(mesh_, u_, degree_) 
        EF_ = mesh_.EFAdj;
        [~, diffs_] = create_local_frames(mesh_);
        r_ = exp(1i*degree_*diffs_);
        E_ = 0;
        for eid_ = 1:mesh_.nE
            if mesh_.isBoundaryEdge(eid_)
                error('Boundary edges not implemented')
            end

            f_ = EF_(eid_, 1);
            g_ = EF_(eid_, 2);
            %A(eid, fi) = r(eid);
            %A(eid, fj) = -1;

            E_ = E_ + abs(u_(f_)*r_(eid_) - u_(g_))^2;
        end
    end

end