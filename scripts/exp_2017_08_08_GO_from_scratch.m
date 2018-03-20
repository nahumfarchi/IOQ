function [] = exp_2017_08_08_GO_from_scratch()
ME = [];
try
%MESHES = {'sphere_s0.off'};
ABOUT = 'Globally optimal from scratch. Following https://github.com/avaxman/DirectionalFieldSynthesis.\r\n\r\n';
OUT_FOLDER_NAME = 'GO_N1_from_scratch';

% MESHES = {'sphere_s0.off', ...
%     'bumpy.off', ...
%     'round_cuber.off', ...
%     'rounded_cube_keenan.off', ...
%     'cow.off', ...
%     'bunny.off', ...
%     'bunny2.off', ...
%     'phands.off', ...
%     'torus_fat_r2.off'};

%MESHES = {'bunny2.off', 'bunny.off', 'rounded_cube_keenan.off', 'round_cuber.off'};
%MESHES = {'rounded_cube_keenan.off'};
%MESHES = {'round_cuber.off'};
MESHES = {'sphere_s0.off'};

PLOT = true;
SAVE = true;
EPS = 1e-10;

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

degree = 1;
soft_ids = [1];


for fname = MESHES
    log_and_print(LOG, 'Loading %s...\r\n', fname{:});
    p = find_data_folder();
    fp = fullfile(p, fname{:});

    mesh = Mesh();
    mesh.loadTM(fp);
    
    V = mesh.V;
    F = mesh.F;
    nF = mesh.nF;
    nE = mesh.nE;
    
    % For convenience, make sure that v0 lies in the same plane 
    % as the face
    %v0 = mesh.V(mesh.F(f0, 2), :) - mesh.V(mesh.F(f0, 1), :);
    f0 = soft_ids(1);
    v0 = V(F(f0, 2), :) - V(F(f0, 1), :);
    soft_vecs = 0.087047269543078509 * v0 / norm(v0);
    
    T1 = zeros(nF, 3);
    T2 = zeros(nF, 3);
    for i = 1:mesh.nF
       e1 = V(F(i,2),:) - V(F(i,1),:);
       e2 = V(F(i,3),:) - V(F(i,2),:);
       T1(i,:) = e1 / norm(e1);
       T2(i,:) = cross(T1(i,:), cross(T1(i,:),e2));
       T2(i,:) = T2(i,:) / norm(T2(i,:));
    end
    
    % Build system matrix
    EF = mesh.EFAdj;
    EV = mesh.EVAdj;
    count = 1;
    I = [];
    J = [];
    S = [];
    frame_diffs = zeros(nE, 1);
    for eid = 1:nE
        f = EF(eid, 1); % Face on the left
        g = EF(eid, 2); % Face on the right
        
        if f <= 0 || g <= 0
            continue
        end
        
        % Compute the complex representation of the common edge
        e = V(EV(eid,2),:) - V(EV(eid,1),:);
        vef = [dot(e, T1(f,:)), dot(e, T2(f,:))];
        vef = vef / norm(vef);
        veg = [dot(e, T1(g,:)), dot(e, T2(g,:))];
        veg = veg / norm(veg);
        cef = vef(1) + 1i*vef(2);
        ceg = veg(1) + 1i*veg(2);
        
        frame_diffs(eid) = angle(ceg) - angle(cef);
        
        % Add the term conj(f)^n*ui - conj(g)^n*uj to the energy matrix
        I(end+1) = count;
        J(end+1) = f;
        S(end+1) = conj(cef)^degree;
        I(end+1) = count;
        J(end+1) = g;
        S(end+1) = -conj(ceg)^degree;
        %S(end+1) = -1;
        
        count = count + 1;
    end
    
    % Convert the constraints into the complex polynomial coefficients and 
    % add them as soft constraints
    lambda = 10e6;
    Ib = [];
    Jb = [];
    Sb = [];
    for i = 1:length(soft_ids)
        f = soft_ids(i);
        v = soft_vecs(i, :);
        c = dot(v, T1(f,:)) + 1i*dot(v, T2(f,:));
        c = 0.03749394817969099-0.078558455845338576*1i;
        
        I(end+1) = count;
        J(end+1) = f;
        S(end+1) = sqrt(lambda);
        Ib(end+1) = count;
        Jb(end+1) = 1;
        Sb(end+1) = c^degree * sqrt(lambda);
        
        count = count + 1;
    end
    
    % Solve with backslash
    A = sparse(I, J, S, count-1, nF);
    b = sparse(Ib, Jb, Sb, count-1, 1);
    %u = A \ b;
    cvx_begin
        variable u(nF) complex
        %minimize ( pow_p(norm(u) - 1, 2) )
        %subject to 
        %    A*u == b
        minimize ( norm(A*u - b) )
    cvx_end
    
    % Convert the interpolated polyvector into Euclidean vectors
    R = zeros(nF, 3);
    for f = 1:nF
        root = find_root(u(f), degree);
        R(f,:) = T1(f,:)*real(root) + T2(f,:)*imag(root);
        R(f,:) = R(f,:) / norm(R(f,:));
    end
    
    GO_bs.ffield = R;
    GO_bs.degree = degree;
    GO_bs.E_MIQ = inf;
    GO_bs.E_GO = norm(A(1:end-1,:)*u-b(1:end-1))^2;
    GO_bs.E_GOc = norm(A*u-b)^2;
    
    % Solve with eig
    nC = length(soft_ids);
    [V, D] = eig(full(A(1:end-nC,:)'*A(1:end-nC,:)));
    u2 = V(:, 1);
    
    % Convert the interpolated polyvector into Euclidean vectors
    R2 = zeros(nF, 3);
    for f = 1:nF
        root = find_root(u2(f), degree);
        R2(f,:) = T1(f,:)*real(root) + T2(f,:)*imag(root);
        R2(f,:) = R2(f,:) / norm(R2(f,:));
    end
    
    GO_eig.u = u2;
    GO_eig.ffield = R2;
    GO_eig.degree = degree;
    GO_eig.E_MIQ = inf;
    GO_eig.E_GO = norm(A(1:end-nC,:)*u2-b(1:end-nC))^2;
    GO_eig.E_GOc = norm(A*u2-b)^2;
    
    % Rotate eig solution so that the constrained face is in the 
    % correct direction
    u3 = u2 * exp(1i * (-angle(u2(1))+angle(u(1))));
    GO_rot.u = u3;
    R3 = zeros(nF, 3);
    for f = 1:nF
        root = find_root(u3(f), degree);
        R3(f,:) = T1(f,:)*real(root) + T2(f,:)*imag(root);
        R3(f,:) = R3(f,:) / norm(R3(f,:));
    end
    
    GO_rot.ffield = R3;
    GO_rot.degree = degree;
    GO_rot.E_MIQ = inf;
    GO_rot.E_GO = norm(A(1:end-1,:)*u3-b(1:end-1))^2;
    GO_rot.E_GOc = norm(A*u3-b)^2;
    
    % Normalize backslash solution
    u4 = u / norm(u);
    GO_nrm.u = u3;
    R4 = zeros(nF, 3);
    for f = 1:nF
        root = find_root(u4(f), degree);
        R4(f,:) = T1(f,:)*real(root) + T2(f,:)*imag(root);
        R4(f,:) = R4(f,:) / norm(R4(f,:));
    end
    
    GO_nrm.ffield = R4;
    GO_nrm.degree = degree;
    GO.nrm.E_MIQ = inf;
    GO_nrm.E_GO = norm(A(1:end-1,:)*u4-b(1:end-1))^2;
    GO_nrm.E_GOc = norm(A*u4-b)^2;
    
    % Globally optimal with svd
    A2 = A(1:end-1, 2:end);
    [U, S, V] = svd(full(A2'*A2));
    u5 = [u(1); V(:, end)];
    GO_SVD.u = u5;
    GO_SVD.E_GO = norm(A(1:end-1,:)*u5-b(1:end-1))^2;
    GO_SVD.E_GOc = norm(A*u5-b)^2;
    R5 = zeros(nF, 3);
    for f = 1:nF
        root = find_root(u5(f), degree);
        R5(f,:) = T1(f,:)*real(root) + T2(f,:)*imag(root);
        R5(f,:) = R5(f,:) / norm(R5(f,:));
    end
    GO_SVD.ffield = R5;
    GO_SVD.degree = degree;
    
    GO_SVD.theta = angle(u5);
    GO_SVD.ffield = angles_to_ffield(GO_SVD.theta, [T1; T2], degree);
    GO_SVD.p = zeros(nE, 1);
    for eid = 1:nE
        fi = EF(eid, 1);
        fj = EF(eid, 2);
        GO_SVD.p(eid) = round(-degree/(2*pi) * (GO_SVD.theta(fi) + frame_diffs(eid) - GO_SVD.theta(fj)));
    end
    GO_SVD.E_MIQ = E_MIQ(mesh, GO_SVD.theta, frame_diffs, GO_SVD.p, degree);
    Kg = get_gaussian_curvature(mesh);
    [d0,d1] = get_exterior_derivatives(mesh);
    [GO_SVD.x, GO_SVD.k] = MIQ_to_TCoDS(GO_SVD.theta, GO_SVD.p, frame_diffs, Kg, d0, d1, degree);
    inds = find(abs(GO_SVD.k) > EPS);
    GO_SVD.S = [inds, GO_SVD.k(inds)];
    
%     GO.theta = angle(u) / degree;
%     %GO.ffield = R; % TODO only works for degree=1 atm
%     GO.ffield = angles_to_ffield(GO.theta, [T1; T2], degree);
%     GO.degree = degree;
%     GO.p = zeros(nE, 1);
%     for eid = 1:nE
%         f = EF(eid, 1);
%         g = EF(eid, 2);
%         GO.p(eid) = round(-degree/(2*pi) * (GO.theta(f) + frame_diffs(eid) - GO.theta(g)));
%     end
%     GO.E = E_MIQ(mesh, GO.theta, frame_diffs, GO.p, degree);
%     Kg = get_gaussian_curvature(mesh);
%     [d0, d1] = get_exterior_derivatives(mesh);
%     [GO.x, GO.k] = MIQ_to_TCoDS(GO.theta, GO.p, frame_diffs, Kg, d0, d1, degree);
%     %inds = find(abs(GO.k) > EPS);
%     %GO.S = [inds, GO.k(inds)];
    
    %GO.ffield = R;
    %GO.degree = degree;
    
    GO_bs.title = {'Backslash', ...
            sprintf('$E_{GO} = %g$', GO_bs.E_GO), ...
            sprintf('$E_{GOc} = %g$', GO_bs.E_GOc)};
    GO_eig.title = {'Eig', ...
            sprintf('$E_{GO} = %g$', GO_eig.E_GO), ...
            sprintf('$E_{GOc} = %g$', GO_eig.E_GOc)};
    GO_rot.title = {'Rotated eig', ...
            sprintf('$E_{GO} = %g$', GO_rot.E_GO), ...
            sprintf('$E_{GOc} = %g$', GO_rot.E_GOc)};
    GO_nrm.title = {'Normalized backslash', ...
            sprintf('$E_{GO} = %g$', GO_nrm.E_GO), ...
            sprintf('$E_{GOc} = %g$', GO_nrm.E_GOc)};
    GO_SVD.title = {'GO SVD', ...
            sprintf('$E_{GO} = %g$', GO_SVD.E_GO), ...
            sprintf('$E_{GOc} = %g$', GO_SVD.E_GOc), ...
            sprintf('$E_{MIQ} = %g$', GO_SVD.E_MIQ)};
    
    if SAVE && PLOT
        MeshVis.wfigs(fname{:}, ...
            mesh, ...
            'OutFolder', OUT_FOLDER, ...
            'Titles', {GO_bs.title, GO_eig.title, GO_rot.title, GO_SVD.title}, ...
            'Nrosy', {GO_bs, GO_eig, GO_rot, GO_SVD});
    elseif PLOT
        figure(2)
        subplot(221);
        MeshVis.plot(mesh, ...
            'nrosy', GO_bs, ...
            'Title', GO_bs.title);
        subplot(222);
        MeshVis.plot(mesh, ...
            'nrosy', GO_eig, ...
            'Title', GO_eig.title);
        subplot(223);
        MeshVis.plot(mesh, ...
            'nrosy', GO_rot, ...
            'Title', GO_rot.title);
        subplot(224);
        MeshVis.plot(mesh, ...
            'nrosy', GO_nrm, ...
            'Title', GO_nrm.title);
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