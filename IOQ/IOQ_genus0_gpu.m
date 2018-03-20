function [k,Lp,E_hist,m_hist] = IOQ_cpu_genus0(V, F, varargin)
	% Input:
	%	V - vertices (nv x 3)
	%	F - faces (nf x 3)
	%
	% Optional name-value pairs:
	%	'Iterations' <n> 
	%		Maximum number of iterations (1000 by default).
	%	'NSingularities' <n>
	%		Number of starting singularities (4*xi by default).
	%	'Laplacian' <'conn' | 'cot' | L>
	%		Laplacian type ('conn' by default).
	%	'LaplacianPInv' <Lp>
	%		The pseudo-inv of the Laplacian. The Laplacian argument is ignored
	%		if this is given since L is only needed to compute Lp.
	%	'Plot' <true | false>
	%		Plots the energy history and m (false by default).
	%	'Tol' <x>
	%		Stopping tolerance (1e-10 by default)
	%
	% Output:
	%	k - singularities (integer vector)

	%%%%% Parse input

	p = inputParser;

	addOptional(p, 'Iterations', 1000);
    addOptional(p, 'NSingularities', []);
    addOptional(p, 'Laplacian', 'conn');
    addOptional(p, 'LaplacianPInv', []);
    addOptional(p, 'Plot', false);
    addOptional(p, 'Tol', 1e-7);
    addOptional(p, 'BlockSize', 8000);

    parse(p, varargin{:});
    opt = p.Results;
    mesh = Mesh(V, F);
    nv = mesh.nV;

    %%%%% Setup

    if isempty(opt.LaplacianPInv)
    	if ischar(opt.Laplacian)
    		disp('Creating L...')
    		if strcmpi(opt.Laplacian, 'conn')
    			[d0, ~] = get_exterior_derivatives(mesh);
    			L = d0'*d0;
    		elseif strcmpi(opt.Laplacian, 'cot')
    			L = -cotmatrix(mesh.V, mesh.F);
    		end
    	else
    		L = opt.Laplacian;
    	end

    	disp('Inverting L...')
    	%Lp = invChol_mex(full(L + 1/nv)) - 1/nv;
        Lp = block_inv_gpu(full(L+1/nv), opt.BlockSize) - 1/nv;
    else
    	Lp = opt.LaplacianPInv;
    end

    disp('Setting initial singularities...')
    Kg = get_gaussian_curvature(mesh);
    k0 = (2 / pi) * Kg;
    if isempty(opt.NSingularities)
    	c = round(abs(sum(k0))); % = 4 xi
    else
    	c = opt.NSingularities;
    end

    k = zeros(nv, 1);
    n_pos_sing = (c + round(sum(k0))) / 2;
    n_neg_sing = (c - round(sum(k0))) / 2;
    inds_pos = randperm(nv, n_pos_sing);
    inds_neg = randperm(nv, n_neg_sing);
    k(inds_pos) = 1;
    k(inds_neg) = -1;

    dlp = diag(Lp);
    b = 2*(Lp*(k - k0));

    assert(n_pos_sing + n_neg_sing == c)
    assert(abs(sum(k)-sum(k0)) < 1e-10)
    
    E_hist = [];
    m_hist = [];
    
    % Setup the gpu kernel
    kernel = parallel.gpu.CUDAKernel('lattice_inner_loop.ptx', 'lattice_inner_loop.cu', 'reduce_cols');
    num_elem = nv^2;
    kernel.ThreadBlockSize = [1, kernel.MaxThreadsPerBlock, 1];
    kernel.GridSize = [1, ceil(nv / kernel.MaxThreadsPerBlock)];
    out_min = single(zeros(1, nv, 'gpuArray'));
    out_r = uint32(zeros(1, nv, 'gpuArray'));
    Lp = gpuArray(single(Lp));
    dlp = gpuArray(single(dlp));
    b = gpuArray(single(b));
    k = gpuArray(single(k));
    %k0 = gpuArray(single(k0));

    %%%%% Minimize energy
    disp('Minimizing...')
    
    for iter = 1:opt.Iterations       
    	%[m, i] = min(bsxfun(@plus, dlp+b, (dlp-b)') - 2*Lp);
        %[m, j] = min(m);
        %i = i(j);
        disp(iter)
        [out_min, out_r] = feval(kernel, out_min, out_r, Lp, dlp, b, nv, nv, num_elem);
        out_r = out_r + 1;
        [m, j] = min(out_min);
        i = out_r(j);

        if abs(m(1)) < opt.Tol
        	disp('Minimization complete because optimality tolerance was reached.')
        	break
        end

        k(i) = k(i) + 1;
        k(j) = k(j) - 1;
        b = b + 2*Lp(:, i) - 2*Lp(:,j);

        if opt.Plot
        	E_hist(iter) = gather((k-k0)'*Lp*(k-k0));
        	m_hist(iter) = gather(m);
        end
    end
    
    if iter == opt.Iterations
        fprintf('Optimization terminated because the iterations budget %d was reached\n', opt.Iterations);
        fprintf('abs(m) = %g\n', abs(m));
    end
    
    k = double(gather(k));
    if nargout > 1
        Lp = double(gather(Lp));
    end

    assert(abs(sum(k)-sum(k0)) < 1e-7)

    if opt.Plot
    	disp('Plotting...')
    	xx = 1:length(E_hist);
    	figure()
    	plot(xx, E_hist, 'r', ...
    		 xx, m_hist, 'b')
    	legend('Energy', 'm')
    	xlabel('Iteration')
    	title('IOQ (cpu)')
    end

    fprintf('Done.\n\n');
end