function [GO_energy, MIQ_energy, x, ffield, S, elapsed] = ...
    calc_GO_and_MIQ_energies(...
        mesh_filename, ...
        degree, ...
        face_ids, constraint_vectors, ...
        eps)

    % Compute GODF and MIQ energy of the given mesh.
    %
    % Input:
    %   mesh_filename - path to .off file
    %   degree - field degree. Default is 1.
    %   face_ids - ids of constrained faces. Default is [1].
    %   constraint_vectors - nCx3 matrix of constraint vectors, where
    %       nC is the number of constraints. Default is the first edge of the face.
    %   eps - will use lambda_1-eps in GODF system. Default is 1e-10.
    %
    % Output:
    %   GO_energy
    %   MIQ_energy
    %   x - solution to (A_s-(lambda1-eps)*I)*x = phi0,
    %       where lambda_1 is the first eigenvalue of A_s.
    %       x(1:nF) - real components of the field
    %       x(nF+1:end) - imaginary components of the field
    %   ffield - (nF*degree)x3 matrix of vector field on faces. First nF
    %       rows are the first vector on each face, second nF rows are the
    %       second vector, and so on until nF*degree.
    %   S - nSx2 matrix of singularity pairs (fid, ki), where fid is the
    %       face id and ki is the singularity index.
	%
	% Examples:
	%	[GO_energy, MIQ_energy, x] = calc_GO_and_MIQ_energies('../data/sphere_s0.off', 4);
	%	[GO_energy, MIQ_energy, x] = calc_GO_and_MIQ_energies('../data/sphere_s0.off', 4, [1,2], [1,0,0;0,1,0], 1e-5);
    
    %%%%%%%%%
    % Setup %
    %%%%%%%%%
    if nargin < 2
        degree = 1;
    end
    
    mesh = Mesh();
    mesh.loadTM(mesh_filename);
    nF = mesh.nF;
    nE = mesh.nE;
    EF = mesh.EFAdj;
    Kg = get_gaussian_curvature(mesh);
    [d0, d1] = get_exterior_derivatives(mesh);
    [local_frames, diffs] = create_local_frames(mesh);
    
    tic
    
    if nargin < 3
        face_ids = [1];
    end
    if nargin < 4
        % Use some arbitrary constraints
        constraint_vectors = zeros(length(face_ids), 3);
        for k = 1:length(face_ids)
            fid = face_ids(k);
            constraint_vectors(k, :) = mesh.V(mesh.F(fid, 2), :) ...
                - mesh.V(mesh.F(fid, 1), :);
        end
        constraint_vectors = normalize_rows(constraint_vectors);
    end
    if nargin < 5
        eps = 1e-10;
    end
       
    %%%%%%%%%%%%%%%%%%%%
    % Globally optimal %
    %%%%%%%%%%%%%%%%%%%%
    [A_s, A_l] = build_real_GO_system(mesh, face_ids, constraint_vectors, degree);
    phi0 = A_l(:);
    
    [V,D] = eigs(A_s, 6, 'sm');
    %lamt = D(1,1) - eps;

    %x = (A_s - lamt*speye(size(A_s))) \ phi0;
    x = V(:, 1);
    x = x / norm(x);
    
    elapsed = toc;
    
    GO_energy = full(x'*A_s*x);
    MIQ_energy = calc_MIQ_energy(mesh, x, degree);
    
    %DD = diag(D) - D(1,1);
    %n_lam1 = length(find(abs(DD)<EPS));
    %U = V(:, 1:n_lam1);
    %z = U*(U'*x);
    
    u = x(1:nF) + 1i*x(nF+1:end);
    theta = angle(u) / degree;
    ffield = angles_to_ffield(theta, local_frames, degree);
    
    % Find singularities
    p = zeros(mesh.nE, 1);
    for eid = 1:nE
        fi = EF(eid, 1);
        fj = EF(eid, 2);
        p(eid) = round(-degree/(2*pi) * (theta(fi) + diffs(eid) - theta(fj)));
    end
    [~, k] = MIQ_to_TCoDS(theta, p, diffs, Kg, d0, d1, degree);
    inds = find(abs(k) > 1e-10);
    S = [inds, k(inds)];

function [A_s, A_l] = build_real_GO_system(mesh, face_ids, constraint_vectors, degree)
    EF = mesh.EFAdj;
    nF = mesh.nF;
    [local_frames, frame_diffs] = create_local_frames(mesh);
    
    % Construct smooth energy quadratic form 
    %   E_s = sum |r_{fg}*u_f - u_g|^2
    %       = x'*A_s*x
    % where x is the real representation of the field, and u is the complex
    % representation:
    %   u_f = x(f)+1i*x(nF+f)
    %   u_g = x(g)+1i*x(nF+g)
    I = [];
    J = [];
    S = [];
    count = 1;
    for eid = 1:mesh.nE
        if mesh.isBoundaryEdge(eid)
            error('Boundary edges not implemented')
        end
        
        f = EF(eid, 1); % Face to the left
        g = EF(eid, 2); % Face to the right
        
        % See notes/2017-08-02 globally optimal.pdf
        I(end+1) = count;
        J(end+1) = f;
        S(end+1) = cos(degree*frame_diffs(eid));
        
        I(end+1) = count;
        J(end+1) = nF+f;
        S(end+1) = -sin(degree*frame_diffs(eid));
        
        I(end+1) = count;
        J(end+1) = g;
        S(end+1) = -1;
        
        I(end+1) = count+1;
        J(end+1) = f;
        S(end+1) = sin(degree*frame_diffs(eid));
        
        I(end+1) = count+1;
        J(end+1) = nF+f;
        S(end+1) = cos(degree*frame_diffs(eid));
        
        I(end+1) = count+1;
        J(end+1) = nF+g;
        S(end+1) = -1;
        
        count = count + 2;
    end
    
    % Constraint energy E_l(psi) = -<psi,phi0>
    k = 1;
    for fid = face_ids
        v = constraint_vectors(k, :);
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
    
    assert(2*mesh.nE + length(face_ids) == count-1);
    n = 2*mesh.nE + length(face_ids);
    m = 2*mesh.nF;
    X = sparse(I, J, S, n, m);
    
    A_s = X(1:end-1,:)'*X(1:end-1,:);
    A_l = X(end,:);
    
function E = calc_MIQ_energy(mesh, x, degree)
    EF = mesh.EFAdj;
    nE = mesh.nE;
    nF = mesh.nF;
    [~, frame_diffs] = create_local_frames(mesh);
    
    u = x(1:nF) + 1i*x(nF+1:end);
    theta = angle(u) / degree;
    periods = zeros(mesh.nE, 1);
    for eid = 1:nE
        fi = EF(eid, 1); 
        fj = EF(eid, 2); 
        periods(eid) = round(-degree/(2*pi) * ...
            (theta(fi) + frame_diffs(eid) - theta(fj)));
    end
    
    E = E_MIQ(mesh, theta, frame_diffs, periods, degree);

