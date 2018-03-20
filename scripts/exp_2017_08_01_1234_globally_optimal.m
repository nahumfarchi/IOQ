function [] = exp_2017_08_01_1234_globally_optimal()
ME = [];
try
%MESHES = {'sphere_s0.off'};
ABOUT = 'Globally optimal, backslash\r\n\r\n';
OUT_FOLDER_NAME = 'GO_N1_backslash';

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
MESHES = {'sphere_s0.off'};

%VIEW_ANGLE = [];

VERBOSE = true;
EPS = 1e-9;

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
    LOG = -1;
end

log_and_print(LOG, '%s\r\n', mfilename('fullpath'));
log_and_print(LOG, ABOUT);
log_and_print(LOG, 'logid : %d\r\n', LOG);

SOLVE_WITH_EIG = false;
SOLVE_WITH_BACKSLASH = ~SOLVE_WITH_EIG;

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
    
    %%%%%%
    % TC %
    %%%%%%
    log_and_print(LOG, '\r\nSolving with TCoDS...\r\n');
    TC_MIQ = TCODS(mesh, MIQ.S, f0, theta0, degree, VERBOSE);
    
    %%%%%%%%%%%%%%%%%%%%
    % Globally optimal %
    %%%%%%%%%%%%%%%%%%%%
    log_and_print(LOG, '\r\nSolving with GO...\r\n');
    EF = mesh.EFAdj;
    [local_frames, diffs] = create_local_frames(mesh);
    r = exp(1i*degree*diffs);
    count = 1;
    I = [];
    J = [];
    S = [];
    for eid = 1:mesh.nE
        if mesh.isBoundaryEdge(eid)
            error('Boundary edges not implemented')
        end
        
        f = EF(eid, 1);
        g = EF(eid, 2);
        %A(eid, fi) = r(eid);
        %A(eid, fj) = -1;
        
        I(end+1) = count;
        J(end+1) = f;
        S(end+1) = r(eid);
        I(end+1) = count;
        J(end+1) = g;
        S(end+1) = -1;
        
        count = count + 1;
    end
    
    %nC = length(f0);
    %C = sparse(nC, mesh.nF);
    %C(1, f0) = 1;
    %d = [i];
    
    % Convert the constraints into the complex polynomial coefficients and 
    % add them as soft constraints
    lambda = 10e6;
    Ib = [];
    Jb = [];
    Sb = [];
    v0 = 0.0880 * v0 / norm(v0);
    for i = 1:length(f0)
        if SOLVE_WITH_EIG
            break;
        end
        f = f0(i);
        v = v0(i, :);
        c = dot(v, local_frames(f,:)) + ...
            1i*dot(v, local_frames(f+mesh.nF,:));
        %c = -0.0844 + 1i*0.025;
        
        I(end+1) = count;
        J(end+1) = f;
        S(end+1) = sqrt(lambda);
        Ib(end+1) = count;
        Jb(end+1) = 1;
        Sb(end+1) = c^degree * sqrt(lambda);
        
        count = count + 1;
    end
    
    % Solve with backslash
    %A = [A; sqrt(lambda)*C];
    %b = [zeros(mesh.nE, 1); sqrt(lambda)*d];
    A = sparse(I, J, S, count-1, mesh.nF);
    b = sparse(Ib, Jb, Sb, count-1, 1);
    if SOLVE_WITH_BACKSLASH
        u = A \ b;
    else
        [V, D] = eig(full(A'*A));
        u = V(:, 1);
        check_norm('norm(A*u-b)^2', 'D(1)');
        figure
        plot(diag(D), '.');
        title('Eigenvalues')
    end
    
    GO.u = u;
    GO.E_GO = norm(A*u-b)^2;
    %check_norm('GO.E_GO', 'D(1)');
    
    % Try rotating u by some amount and see if the result still has
    % the same energy and is still an eigenvector of A'*A.
    %u_rot = exp(1i*0.32421241241)*u;
    %check_norm('norm(A*u-b)', 'norm(A*u_rot-b)');
    %check_norm('(A''*A)*u', 'D(1)*u');
    %check_norm('(A''*A)*u_rot', 'D(1)*u_rot');
    
    GO.theta = angle(u) / degree;
    GO.ffield = angles_to_ffield(GO.theta, local_frames, degree);
    GO.degree = degree;
    GO.p = zeros(mesh.nE, 1);
    for eid = 1:mesh.nE
        fi = EF(eid, 1);
        fj = EF(eid, 2);
        GO.p(eid) = round(-degree/(2*pi) * (GO.theta(fi) + diffs(eid) - GO.theta(fj)));
    end
    GO.E_MIQ = E_MIQ(mesh, GO.theta, diffs, GO.p, degree);
    [GO.x, GO.k] = MIQ_to_TCoDS(GO.theta, GO.p, diffs, MIQ.Kg, d0, d1, degree);
    inds = find(abs(GO.k) > EPS);
    GO.S = [inds, GO.k(inds)];
    
    %%%%%%
    % TC %
    %%%%%%
    f0 = [1];
    log_and_print(LOG, '\r\nSolving with TCoDS...\r\n');
    I = find(abs(GO.k) > EPS);
    TC_GO = TCODS(mesh, [I, GO.k(I)/degree], f0, theta0, degree, VERBOSE);
   
    % Calculate GO energy
    MIQ.u = langles_to_complex(MIQ.theta);
%     for fid = 1:mesh.nF
%         v = MIQ.ffield(fid, :);
%         c = dot(v, local_frames(fid,:)) + 1i*dot(v, local_frames(fid+mesh.nF,:));
%         MIQ.u2(fid) = c;
%     end
    MIQ.E_GO = norm(A*MIQ.u-b)^2;
    
    TC_GO.u = langles_to_complex(TC_GO.theta);
    %for fid = 1:mesh.nF
    %    v = TC_GO.ffield(fid, :);
    %  c = dot(v, local_frames(fid,:)) + 1i*dot(v, local_frames(fid+mesh.nF,:));
    %    TC_GO.u(fid) = c;
    %end
    TC_GO.E_GO = norm(A*TC_GO.u-b)^2;
    
    % Globally optimal with SVD
    A2 = A(1:end-1, 2:end);
    [U, S, V] = svd(full(A2'*A2));
    u = [GO.u(1); V(:, end)];
    GO_SVD.u = u;
    GO_SVD.E_GO = norm(A*u - b)^2;
    
    GO_SVD.theta = angle(u) / degree;
    GO_SVD.ffield = angles_to_ffield(GO_SVD.theta, local_frames, degree);
    GO_SVD.degree = degree;
    GO_SVD.p = zeros(mesh.nE, 1);
    for eid = 1:mesh.nE
        fi = EF(eid, 1);
        fj = EF(eid, 2);
        GO_SVD.p(eid) = round(-degree/(2*pi) * (GO_SVD.theta(fi) + diffs(eid) - GO_SVD.theta(fj)));
    end
    GO_SVD.E_MIQ = E_MIQ(mesh, GO_SVD.theta, diffs, GO_SVD.p, degree);
    [GO_SVD.x, GO_SVD.k] = MIQ_to_TCoDS(GO_SVD.theta, GO_SVD.p, diffs, MIQ.Kg, d0, d1, degree);
    inds = find(abs(GO_SVD.k) > EPS);
    GO_SVD.S = [inds, GO_SVD.k(inds)];
       
    %%%%%%%%%
    % Plots %
    %%%%%%%%%
    
    MIQ.title = {'MIQ', sprintf('$E_{MIQ} = %g$', MIQ.E), sprintf('$E_{GO} = %g$', MIQ.E_GO)};
    TC_MIQ.title = {'TC (MIQ)', sprintf('$E_{MIQ} = %g$', TC_MIQ.E)};
    GO.title = {'GO', sprintf('$E_{MIQ} = %g$', GO.E_MIQ), sprintf('$E_{GO} = %g$', GO.E_GO)};
    TC_GO.title = {'TC (GO) ', sprintf('$E_{MIQ} = %g$', TC_GO.E), sprintf('$E_{GO} = %g$', TC_GO.E_GO)};
    GO_SVD.title = {'GO SVD ', sprintf('$E_{MIQ} = %g$', GO_SVD.E_MIQ), sprintf('$E_{GO} = %g$', GO_SVD.E_GO)};
    
    if SAVE && PLOT
        log_and_print(LOG, '\r\nPlotting...\r\n');
        %MeshVis.wfigs('abc', mesh);
        MeshVis.wfigs(fname{:}, ...
            mesh, ...
            'OutFolder', OUT_FOLDER, ...
            'Titles', {MIQ.title, TC_MIQ.title, GO.title, GO_SVD.title}, ...
            'Nrosy', {MIQ, TC_MIQ, GO, GO_SVD}, ...
            'Funcs', {[], [], full(abs(GO.u)), []});
    elseif PLOT
        figure
        subplot(221)
        MeshVis.plot(mesh, 'nrosy', MIQ);
        title(['MIQ ', num2str(MIQ.E)])
        subplot(222)
        MeshVis.plot(mesh, 'nrosy', TC_MIQ);
        title(['TCoDS ', num2str(TC_MIQ.E)])
        subplot(223)
        MeshVis.plot(mesh, 'nrosy', GO);
        title(['GO ', num2str(GO.E)])
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

end