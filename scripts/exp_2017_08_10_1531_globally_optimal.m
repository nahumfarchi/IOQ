function [] = exp_2017_08_01_1234_globally_optimal()
ME = [];
try
%MESHES = {'sphere_s0.off'};
ABOUT = 'Globally optimal, backslash\r\n\r\n';
OUT_FOLDER_NAME = 'GO_svd';

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
    
    A = sparse(I, J, S, count-1, mesh.nF);
    A2 = A(1:end-1, :);
    b = sparse(Ib, Jb, Sb, count-1, 1);
    b2 = b(1:end-1);
    
    % Solve with backslash
    log_and_print(LOG, '\r\nSolving GO with backslash...\r\n');
    GO_bs.u = A \ b;
    GO_bs.E_GO = norm(A2*GO_bs.u-b2)^2;
    GO_bs.E_GOc = norm(A*GO_bs.u-b)^2;
    GO_bs.theta = angle(GO_bs.u) / degree;
    GO_bs.ffield = angles_to_ffield(GO_bs.theta, local_frames, degree);
    GO_bs.degree = degree;
    GO_bs.p = zeros(mesh.nE, 1);
    for eid = 1:mesh.nE
        fi = EF(eid, 1);
        fj = EF(eid, 2);
        GO_bs.p(eid) = round(-degree/(2*pi) * (GO_bs.theta(fi) + diffs(eid) - GO_bs.theta(fj)));
    end
    GO_bs.E_MIQ = E_MIQ(mesh, GO_bs.theta, diffs, GO_bs.p, degree);
    [GO_bs.x, GO_bs.k] = MIQ_to_TCoDS(GO_bs.theta, GO_bs.p, diffs, MIQ.Kg, d0, d1, degree);
    inds = find(abs(GO_bs.k) > EPS);
    GO_bs.S = [inds, GO_bs.k(inds)];
    
    % Solve with eig
    log_and_print(LOG, '\r\nSolving GO with eig...\r\n');
    [V, D] = eig(full(A2'*A2));
    GO_eig.u = V(:, 1);
    GO_eig.E_GO = norm(A2*GO_eig.u-b2)^2;
    GO_eig.E_GOc = norm(A*GO_eig.u-b)^2;
    GO_eig.theta = angle(GO_eig.u) / degree;
    GO_eig.ffield = angles_to_ffield(GO_eig.theta, local_frames, degree);
    GO_eig.degree = degree;
    GO_eig.p = zeros(mesh.nE, 1);
    for eid = 1:mesh.nE
        fi = EF(eid, 1);
        fj = EF(eid, 2);
        GO_eig.p(eid) = round(-degree/(2*pi) * (GO_eig.theta(fi) + diffs(eid) - GO_eig.theta(fj)));
    end
    GO_eig.E_MIQ = E_MIQ(mesh, GO_eig.theta, diffs, GO_eig.p, degree);
    [GO_eig.x, GO_eig.k] = MIQ_to_TCoDS(GO_eig.theta, GO_eig.p, diffs, MIQ.Kg, d0, d1, degree);
    inds = find(abs(GO_eig.k) > EPS);
    GO_eig.S = [inds, GO_eig.k(inds)];
    
    % Calculate GO energy
    MIQ.u = langles_to_complex(MIQ.theta);
    MIQ.E_GO = norm(A2*MIQ.u-b2)^2;
    MIQ.E_GOc = norm(A*MIQ.u-b)^2;
    
    % Globally optimal with SVD
    log_and_print(LOG, '\r\nSolving with SVD...\r\n');
    A3 = A(1:end-1, 2:end);
    [U, S, V] = svd(full(A3'*A3));
    u = [GO_bs.u(1); V(:, end)];
    GO_SVD.u = u;
    GO_SVD.E_GO = norm(A2*u- b2)^2;
    GO_SVD.E_GOc = norm(A*u - b)^2;
    
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
    
    MIQ.title = {'MIQ', ...
        sprintf('$E_{MIQ} = %.4g, E_{GO} = %.4g, E_{GOc} = %.4g$', MIQ.E, MIQ.E_GO, MIQ.E_GOc)};
    GO_bs.title = {'GO (backslash)', ...
        sprintf('$E_{MIQ} = %.4g, E_{GO} = %.4g, E_{GOc} = %.4g$', GO_bs.E_MIQ, GO_bs.E_GO, GO_bs.E_GOc)};
    GO_eig.title = {'GO (eig)', ...
        sprintf('$E_{MIQ} = %g, E_{GO} = %.4g, E_{GOc} = %.4g$', GO_eig.E_MIQ, GO_eig.E_GO, GO_eig.E_GOc)};
    GO_SVD.title = {'GO (SVD)', ...
        sprintf('$E_{MIQ} = %g, E_{GO} = %.4g, E_{GOc} = %.4g$', GO_SVD.E_MIQ, GO_SVD.E_GO, GO_SVD.E_GOc)};
    
    if SAVE && PLOT
        log_and_print(LOG, '\r\nPlotting...\r\n');
        %MeshVis.wfigs('abc', mesh);
        MeshVis.wfigs(fname{:}, ...
            mesh, ...
            'OutFolder', OUT_FOLDER, ...
            'Titles', {MIQ.title, GO_bs.title, GO_eig.title, GO_SVD.title, GO_SVD.title}, ...
            'Nrosy', {MIQ, GO_bs, GO_eig, GO_SVD});
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