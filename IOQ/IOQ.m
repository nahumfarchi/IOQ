function [alpha, beta, connection, stats, out] = IOQ(verts, faces, varargin)
    % function [alpha, beta, stats, Lp] = IOQ(verts, faces, varargin)
	% Input:
	%	verts - (nv x 3)
	%	faces - (nf x 3)
    %
    % Pair arguments:
    %   'Iterations', n             Max number of iterations (default 2000)
    %   'Laplacian', 'conn | cot'   Laplacian type (def. 'conn')
    %   'Lp', Lp                    Laplacian pseudoinverse (def. [])
    %   'GPU', true | false         (def. true)
    %   'Verbose', true | false     (def. true)
    %   'InvMethod', 'GPUInv | BlockGPU | CholMexInv | ApproxResistance'
    %                               (def. GPUInv)
    %   'JLEps', eps                (def. 0.5)
    %   'JLFac', fac                (def. 24)
    %   'Mesh', m                   (def. [])
    %   'BetaMethod', 'round | cvp' (def. 'round')
    %   'BetaMax', x                (def. 20)
    %
    % Output:
    %   alpha      - cone singularitie indices (nv x 1 integer vector)
    %   beta       - generator singularities indices (2*ng x 1 integer vector)
    %   connection - the trivial connection (ne x 1 vector with the angle
    %   adjustments)
    %   stats      - a struct with various statistics (energy, timing, etc)
    %   out        - a strcut with various matrices that might be useful (L, Lp, d0, B, etc)
    
    % ---------------------------------------------------------------------
    % Parse input
    % ---------------------------------------------------------------------
    
    parser = inputParser;

	addOptional(parser, 'Iterations', 2000);
    addOptional(parser, 'NSingularities', 0); % Number of initial singularities
    addOptional(parser, 'Laplacian', 'conn');
    addOptional(parser, 'Lp', []);
    addOptional(parser, 'Tol', 1e-6);
    addOptional(parser, 'BlockSize', 20000);
    addOptional(parser, 'GPU', true);
    addOptional(parser, 'Verbose', true);
    addOptional(parser, 'InvMethod', 'GPUInv');
    addOptional(parser, 'JLEps', 0.5);
    addOptional(parser, 'JLFac', 24);
    addOptional(parser, 'Mesh', []);
    addOptional(parser, 'KernelEntry', 'reduce_cols');
 	addOptional(parser, 'alpha', []);
    addOptional(parser, 'BetaMethod', 'round');
    addOptional(parser, 'BetaMax', 20);
    addOptional(parser, 'Constraints', []);
    addOptional(parser, 'Gamma', []);
    addOptional(parser, 'GatherStats', true);

    parse(parser, varargin{:});
    opt = parser.Results;
    
    verb = opt.Verbose;
    approx = strcmpi(opt.InvMethod, 'ApproxResistance');
    if approx
        opt.KernelEntry = 'reduce_cols_R_symmetric';
    end
    
    if isempty(opt.Mesh)
        mesh = Mesh(verts, faces);
    else
        mesh = opt.Mesh;
    end
    
    nv = mesh.nV; ne = mesh.nE;
    
    gather_stats = opt.GatherStats;
    if nargout > 3 && gather_stats
        stats.T_setup = 0;
        stats.T_lap = 0;
        stats.T_inv = 0;
        stats.T_gen = 0;
        stats.T_gridsearch = 0;
        stats.T_beta = 0;
        stats.T_total = 0;
        stats.E_hist = [];
        stats.Ea_hist = [];
        stats.Eb_hist = [];
    end
    
    % ---------------------------------------------------------------------
    % Setup
    % ---------------------------------------------------------------------
    
    tic
    if opt.GPU
        if verb, disp('Setting up CUDA kernel...'); end
        kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', opt.KernelEntry);
        num_elem = nv^2;
        kernel.ThreadBlockSize = [1, kernel.MaxThreadsPerBlock, 1];
        kernel.GridSize = [1, ceil(nv / kernel.MaxThreadsPerBlock)];
        out_min = single(zeros(1, nv, 'gpuArray'));
        out_r = uint32(zeros(1, nv, 'gpuArray'));
        gd = gpuDevice();
    else
        if verb, disp('Not using GPU...'); end
        opt.InvMethod = 'CholMexInv';
        gd = [];
    end
        
    if verb, disp('Placing initial singularities...'); end
    alpha_g = mesh.gaussian_cur;
    genus = round(1 - sum(alpha_g) / (4*pi));
    xi = 2 - 2*genus;
    x0 = (2/pi)*alpha_g;
    if isempty(opt.alpha)
        %ns = round(abs(sum((x0)))); % = 4 xi  
        ns = abs(4 * xi);
		alpha = zeros(nv, 1);
		n_pos_sing = (ns + round(sum(x0))) / 2 + opt.NSingularities;
		n_neg_sing = (ns - round(sum(x0))) / 2 + opt.NSingularities;
        inds = randperm(nv, n_pos_sing+n_neg_sing);
        inds_pos = inds(1:n_pos_sing);
        inds_neg = inds(n_pos_sing+1:end);
		alpha(inds_pos) = 1;
		alpha(inds_neg) = -1;
		assert(n_pos_sing + n_neg_sing - 2*opt.NSingularities == ns)
        
        if gather_stats
            stats.init_alpha = alpha;
            stats.n_init_sing = length(inds);
        end
    else
		alpha = opt.alpha;
		assert(length(alpha) == nv)
    end
    if ~isempty(gd), wait(gd); end; stats.T_setup = toc;
    
    % Create Laplacian and it's inverse (or the approximate resistance distance)
    if verb, disp('Creating L...'); end
    tic
    if isempty(opt.Constraints)
        d0 = get_exterior_derivatives(mesh);
    else
        [d0, d1] = get_exterior_derivatives(mesh);
    end
    L = d0'*d0;
    FL = factorize(L);
    if ~isempty(gd), wait(gd); end; stats.T_lap = toc;
    
    if isempty(opt.Lp)
        if verb, disp('Inverting L...'); end
        tic
        switch opt.InvMethod
            case 'GPUInv'
                Lp = inv(gpuArray(single(L+1/nv))) - 1/nv;
            case 'BlockGPU'
                Lp = block_inv_gpu(full(L+1/nv), opt.BlockSize, true) - 1/nv;
            case 'ApproxResistance'
                [Ztilde, Rtilde] = resistance_distance();
                approx = true;
            case 'CholMexInv'
                Lp = invChol_mex(full(L + 1/nv)) - 1/nv;
            otherwise
                error('Unkown inv method')
        end
        if ~isempty(gd), wait(gd); end; stats.T_inv = toc;
    else
        Lp = opt.Lp;
        stats.T_lap = nan;
        stats.T_inv = nan;
    end
    
    % Generator cycles
    tic
    if genus == 0
        beta = [];
        beta_g = [];
        B = [];
    else
        if verb, disp('Generator cycles...'); end
        H = mesh.H;
        %beta_g = wrapToPi(generator_angle_defects(mesh));
        beta_g = mesh.generator_defects;
        B = H - d0 * (FL \ (d0' * H));
        FHB = factorize(H'*B);
        A2 = (H'*d0) / FL;
        y0 = (2/pi) * (A2*alpha_g - beta_g);
    end
    if ~isempty(gd), wait(gd); end; stats.T_gen = toc;
    
    if isempty(opt.Constraints)
        has_constraints = false;
    else
        has_constraints = true;
        constrained_faces = opt.Constraints(:, 1);
        constraint_thetas = opt.Constraints(:, 2);
        nc = length(constrained_faces);
        f0 = opt.Constraints(1, 1); theta0 = opt.Constraints(1, 2);
        %frame1 = local_frames(f0, :); frame2 = local_frames(f0+nf);
        %gvec = cos(theta0)*frame1 + sin(theta0)*frame2;
        [local_frames, frame_diffs] = create_local_frames(mesh);
        [R, gamma_g] = create_constraints_mat(mesh, constrained_faces, constraint_thetas, frame_diffs);
        R = R';
    end
    
    a = []; b = []; c = [];
    
    % ---------------------------------------------------------------------
    % Solve for alpha (gridsearch)
    % ---------------------------------------------------------------------
    
    if verb, disp('Gridsearch...'); end
    tic
    if approx
    elseif has_constraints
        alpha = gridsearch2(alpha);
    else
        alpha = gridsearch(Lp, alpha, x0);
    end
    if ~isempty(gd), wait(gd); end; stats.T_gridsearch = toc;
  
    a = FL \ ( (pi/2)*alpha - alpha_g );
    connection = d0*a;
    
    % ---------------------------------------------------------------------
    % Solve for beta (round or cvp)
    % ---------------------------------------------------------------------
    
    tic
    if genus > 0
        switch opt.BetaMethod
            case 'round'
                beta = round(full(A2)*alpha - y0);
            case 'cvp'
                G = (B / FHB) * (pi/2);
                m = size(G, 2);
                target = B * ( FHB \ ( beta_g + H'*d0*a ) );
                beta = SEA_det(m, opt.BetaMax, target, G, false);
            otherwise
                error('Unkown BetaMethod')
        end
        
        b = FHB \ ( (pi/2)*beta - beta_g - H'*d0*a );
        connection = connection + B*b;
    end
    if ~isempty(gd), wait(gd); end; stats.T_beta = toc;
    
    % ---------------------------------------------------------------------
    % Solve for gamma (i.e., the constraint singularities)
    % ---------------------------------------------------------------------
    if ~isempty(opt.Constraints)
        if ~isempty(opt.Gamma)
            gamma = opt.Gamma;
        elseif genus > 0
            gamma = round( 2/pi * ( gamma_g + R'*d0*a + R'*B*b ) );
        else
            gamma = round( 2/pi * ( gamma_g + R'*d0*a ) );
        end
        
        if genus > 0
            %c = inverse( R' * d1' ) * (pi/2 * gamma - gamma_g - R'*d0*a - R'*B*b );
            c = lsqlin(d1', ...
                zeros(ne, 1), [], [], R'*d1', (pi/2) * gamma - gamma_g - R'*d0*a - R'*B*b, [], []);
        else
            %gamma = zeros(nc - 1, 1);
            %c = (R' * d1') \ (pi/2 * gamma - gamma_g - R'*d0*a );
            c = lsqlin(d1', ...
                zeros(ne, 1), [], [], R'*d1', (pi/2) * gamma - gamma_g - R'*d0*a, [], []);
        end
        connection = connection + d1'*c;
    end
    
    % ---------------------------------------------------------------------
    % Output
    % ---------------------------------------------------------------------
    
    if nargout > 3
        stats.miq_energy = norm(connection)^2;
        stats.T_total = stats.T_setup + ...
            stats.T_lap + ...
            stats.T_inv + ...
            stats.T_gen + ...
            stats.T_gridsearch + ...
            stats.T_beta;
    end
    
    if nargout > 4
        out.L = L;
        out.Lp = gather(Lp);
        out.d0 = d0;
        out.F = FL;
        out.B = B;
        out.alpha_g = alpha_g;
        out.beta_g = beta_g;
        out.gamma = gamma;
        out.a = a;
        out.b = b;
        out.c = c;
        degree = 4;
        [out.ffield, out.theta, local_frames, frame_diffs] = connection_to_nrosy(...
                     mesh, ...
                     connection, ...
                     f0, ...
                     theta0, ...
                     degree, ...
                     'LocalFrames', local_frames, ...
                     'FrameDiffs', frame_diffs);
        out.m = mesh;
        inds = find(alpha);
        vert_sing = [inds, alpha(inds)];
        inds = find(beta);
        gen_sing = [inds, beta(inds)];
        mesh.set_ffield(...
                  degree, ...
                  local_frames, ...
                  frame_diffs, ...
                  out.theta, ...
                  out.ffield, ...
                  vert_sing, ...
                  gen_sing, ...
                  norm(connection)^2, ...
                  connection);
    end
    
    alpha = double(gather(alpha));
    beta  = double(gather(beta));
        
% -------------------------------------------------------------------------
% Helper functions    
% -------------------------------------------------------------------------

function x = gridsearch(M, x, x0)
    % M  - Lp
    % x  - alpha
    % x0 - (2 / pi) * alpha_g 
    if ~isa(M, 'single')
        M = single(M);
    end
    if ~isa(x, 'single')
        x = single(x);
    end
    if ~isa(x0, 'single')
        x0 = single(x0);
    end
    min_val = -inf;
    
    if opt.GPU
        try
            % Calculate b and dm on the cpu first (uses less gpu memory)
            b_ = gpuArray(2*(M*(x - x0)));
            dm = gpuArray(diag(M));
            M = gpuArray(M);
        catch
            % Probably out of memory, so try to pass only the lower
            % triangular part of M, which is symmetric.
            % Note that we are now passing the resistance to the kernel
            % (instead of M).
            if verb, disp('Out of memory on GPU. Trying to use lower triangle part of the resistance distance'); end
            if ~exist('b', 'var')
               b_ = 2*(M*(x - x0));
            end
            if ~exist('dm', 'var')
               dm = diag(M);
            else
               dm = gather(dm);
            end
            % Calculate resistance distance
            Mtril = bsxfun(@plus, dm, dm') - 2*M; 

            % Take just the lower part (not including the diagonal)
            mask = tril(true(size(Mtril)), -1);
            Mtril = Mtril(mask);

            Mtril = gpuArray(single(Mtril));
            dm = gpuArray(dm);

            % Setup the symmetric kernel
            opt.KernelEntry = 'reduce_cols_R_symmetric';
            kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', opt.KernelEntry);
            num_elem = nv^2;
            kernel.ThreadBlockSize = [1, kernel.MaxThreadsPerBlock, 1];
            kernel.GridSize = [1, ceil(nv / kernel.MaxThreadsPerBlock)];
        end
    else
        b_ = 2*(M*(x - x0));
        dm = diag(M);
    end
    
    for iter = 1:opt.Iterations
        if gather_stats
            update_stats(iter, x, x0);
        end

        if opt.Verbose && mod(iter, 10) == 0
            fprintf('iter, m : %d, %g\n', iter, min_val);
        end
        
        if ~opt.GPU || strcmpi(opt.KernelEntry, 'bsx')
            [min_val, i] = min(bsxfun(@plus, dm+b_, (dm-b_)') - 2*M);
            [min_val, j] = min(min_val);
            i = i(j);
        else
            switch opt.KernelEntry
                case 'reduce_cols'
                    [out_min, out_r] = feval(kernel, out_min, out_r, M, dm, b_, nv, nv, num_elem);
                case 'reduce_cols_R_symmetric'
                    [out_min, out_r] = feval(kernel, out_min, out_r, Mtril, b_, nv, nv, num_elem);
                otherwise
                    errror('Unkown kernel entry point')
            end
            [min_val, j] = min(out_min);
            i = out_r(j) + 1;
        end 
        
        if abs(min_val(1)) < opt.Tol
            disp('Minimization complete because optimality tolerance was reached.')
            disp(['m : ', num2str(min_val)])
            break
        end
        
        x(i) = x(i) + 1;
        x(j) = x(j) - 1;
        b_ = b_ + 2*M(:,i) - 2*M(:,j);
    end
    
    if iter == opt.Iterations
        fprintf('Optimization terminated because the iterations budget %d was reached\n', opt.Iterations);
        fprintf('min_val = %g\n', min_val);
    end
end

function x = gridsearch2(alpha)
    x = alpha;
    x0 = (2/pi)*alpha_g;
    M1 = Lp;
    b1 = 2*M1*(x-x0);
    dm1 = diag(M1);
    
    M2 = inverse(R'*d1')'*(d1*d1')*inverse(R'*d1');
    A2 = R'*d0*Lp;
    M2_tilde = A2'*M2*A2;
    dm2 = diag(M2_tilde);

    a = FL \ ( (pi/2)*alpha - alpha_g );
    if ~exist('gamma', 'var')
        gamma = round( 2/pi * ( gamma_g + R'*d0*a ) );
    end
    %y_tilde = gamma - (2/pi)*gamma_g+R'*d0*Lp*x0;
    y_tilde = gamma - (2/pi)*gamma_g+A2*x0;
    b2 = 2*A2'*(M2'*(A2*x-y_tilde));   
    
    for iter = 1:opt.Iterations
        if gather_stats
            update_stats(iter, x, x0);
        end
        
        if opt.Verbose && mod(iter, 10) == 0
            fprintf('iter, m : %d, %g\n', iter, min_val);
        end
        
        [min_val, i] = min(bsxfun(@plus, dm1+b1, (dm1-b1)' - 2*M1) + ...
                           bsxfun(@plus, dm2+b2, (dm2-b2)' - 2*M2_tilde));
        [min_val, j] = min(min_val);
        i = i(j);
        x(i) = x(i) + 1;
        x(j) = x(j) - 1;
        b1 = b1 + 2*M1(:, i) - 2*M1(:, j);
        b2 = b2 + 2*M2_tilde(:, i) - 2*M2_tilde(:, j);
        
%         a = FL \ ( (pi/2)*x - alpha_g );
%         gamma = round( 2/pi * ( gamma_g + R'*d0*a ) );
%         %y_tilde = gamma - (2/pi)*gamma_g+R'*d0*Lp*x0;
%         y_tilde = gamma - (2/pi)*gamma_g+A2*x0;
%         b2 = 2*A2'*(M2'*(A2*x-y_tilde));   
        
        if abs(min_val(1)) < opt.Tol
            disp('Minimization complete because optimality tolerance was reached.')
            disp(['m : ', num2str(min_val)])
            break
        end
    end
end
    
function [Ztilde, Rtilde] = resistance_distance()
%[Ztilde, Rtilde] = resistance_distance(m, eps, JLFac, useGPU)
    eps = opt.JLEps;
    JLFac = opt.JLFac;
    useGPU = opt.UseGPU;
    
    %nv = m.nV; ne = m.nE;
    if ~exist('d0', 'var')
        d0 = get_exterior_derivatives(mesh);
        
    end
    if ~exist('L', 'var')
        L = d0' * d0;
    end
    if ~exist('F', 'var')
        FL = factorize(L);
    end
    %L = d0' * d0;
    %F = factorize(L);
    
    k = round(JLFac * log(nv) / eps^2);
    Y = 2*(rand(k, ne) > 0.5) - 1;
    Y = (1/sqrt(k)) * Y;
    %Y = randi([-1, 1], k, ne) / sqrt(k); % slightly faster
    Y = Y * d0;
    Ztilde = (FL \ Y');
    % note that Rtilde is only the upper triangle part (1d vector)
    if useGPU
        Ztilde = gpuArray(single(Ztilde));
        %Ztilde = gpuArray(Ztilde);
        if nv < 60000
            Rtilde = pdist(Ztilde, 'squaredeuclidean')'; 
        else
            Rtilde = pdist_block(Ztilde);
        end
        %Ztilde = gather(Ztilde); % wasteful
    else
        Rtilde = pdist(Ztilde, 'squaredeuclidean')';
    end

    %Rtilde2 = pdist_block();
    %assert(norm(squareform(Rtilde) - Rtilde2, 'fro') < 1e-4)
    %assert(norm(Rtilde - Rtilde2, 'fro') < 1e-10)
end

    
function res = pdist_block(Ztilde)
    %inds_cell = mat2tiles(1:nv, 1, block_size);
    %Rtilde = zeros(nv, nv);
    %for i = 1:numel(inds_cell)
    %    inds = inds_cell{i};
    %    X = pdist2(Ztilde, Ztilde(inds, :), 'squaredeuclidean');
    %    Rtilde(:, inds) = gather(X);
    %end
    %Rtilde = tril(Rtilde, -1);
    %Rtilde = Rtilde(Rtilde>0);
    %Rtilde = Rtilde(:);
    
    sz = ceil(nv / 2);
    D11 = pdist2(Ztilde(1:sz, :), Ztilde(1:sz, :), 'squaredeuclidean');
    D11 = gather(D11);
    D12 = pdist2(Ztilde(1:sz, :), Ztilde(sz+1:end, :), 'squaredeuclidean');
    D12 = gather(D12);
    D22 = pdist2(Ztilde(sz+1:end, :), Ztilde(sz+1:end, :), 'squaredeuclidean');
    D22 = gather(D22);
    %Rtilde = [D11, D12; D12', D22];
    %Rtilde = tril(Rtilde, -1);
    %Rtilde = Rtilde(Rtilde>0);
    %Rtilde = Rtilde(:);
    x = [D11; D12'];
    mask = tril(true(size(x)), -1);
    x = x(mask);
    mask = tril(true(size(D22)), -1);
    y = D22(mask);
    res = [x(:); y(:)];
    res = gpuArray(res);
end

function update_stats(iter, x, x0)
    %a = FL \ gather(x - x0);
    a = FL \ ( (pi/2)*x - alpha_g );
    stats.Ea_hist(iter) = norm(d0*a)^2;
    
    if genus > 0
        %c = inverse( R' * d1' ) * (pi/2 * gamma - gamma_g - R'*d0*a - R'*B*b );
        [c,~,~,~,~] = lsqlin(d1', ...
            zeros(ne, 1), [], [], R'*d1', (pi/2) * gamma - gamma_g - R'*d0*a - R'*B*b, [], []);
    else
        %gamma = zeros(nc - 1, 1);
        %c = (R' * d1') \ (pi/2 * gamma - gamma_g - R'*d0*a );
        [c,~,~,~,~] = lsqlin(d1', ...
            zeros(ne, 1), [], [], R'*d1', (pi/2) * gamma - gamma_g - R'*d0*a, [], []);
    end
    stats.Ec_hist(iter) = norm(d1'*c)^2;
    
    stats.Ehist(iter) = stats.Ea_hist(iter) + stats.Ec_hist(iter);
    
    %stats.Ea_hist(iter) = gather((x-x0)' * a);
    %if exist('c', 'var')
    %    stats.Ec_hist(iter) = gather(
    % TODO
%     if genus > 0        
%         stats.Eb_hist(iter) = nan;
%     else
%         stats.Eb_hist(iter) = 0;
%     end
%     stats.E_hist(iter) = stats.Ea_hist(iter) + stats.Eb_hist(iter);
end

end
    