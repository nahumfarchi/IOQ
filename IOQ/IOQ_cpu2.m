function [k,Lp,E_hist] = IOQ_cpu(V, F, varargin)
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
	%		The pseudo-inv of the Laplacian. The Laplacian argument value is ignored
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
    addOptional(p, 'Tol', 1e-10);

    parse(p, varargin{:});
    opt = p.Results;
    mesh = Mesh(V, F);
    nv = mesh.nV;

    %%%%% Setup
    [~, Kf, d0, ~, H] = tcods_gsystem(V, F);
    P_H = speye(size(H,1)) - H*inv(H'*H)*H'; 

    if isempty(opt.LaplacianPInv)
    	if ischar(opt.Laplacian)
    		disp('Creating L...')
    		if strcmp(lower(opt.Laplacian), 'conn')
    			[d0, ~] = get_exterior_derivatives(mesh);
    			L = d0'*d0;
            elseif strcmp(lower(opt.Laplacian), 'conn2')
    			L = d0'*P_H*d0;
    		elseif strcmp(lower(opt.Laplacian), 'cot')
    			L = -cotmatrix(mesh.V, mesh.F);
    		end
    	else
    		L = opt.Laplacian;
    	end

    	disp('Inverting L...')
    	Lp = invChol_mex(full(L + 1/nv)) - 1/nv;
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

%     B = H - d0*(d0\H);
%     M = B*inv(H'*B);
%     b_h_0 = Kf(nv+1:end)/(pi/2) - H'*d0*Lp*Kf(1:nv)/(pi/2);
%     b_h = b_h_0 + H'*d0*Lp*k;
%     Hd0L = H'*d0*Lp;
%     MM = M'*M;
    
    %%%%%%%%%%%%
    %%% Try to use pinv(d0'*P_H*d0) instead of pinv(d0'*d0)
    %%%%% Minimize energy

    disp('Minimizing...')
    E_hist = [];
    for iter = 1:opt.Iterations
    	[m, i] = min(bsxfun(@plus, dlp+b, (dlp-b)') - 2*Lp);% + ...
%            dot(b_h + Hd0L,MM*(b_h + Hd0L),1)' + dot(Hd0L,MM*Hd0L,1) - 2*(b_h + Hd0L)'*MM*Hd0L - b_h'*MM*b_h);
        [m, j] = min(m);
        i = i(j);

%        assert(m < 0);
        if abs(m(1)) < opt.Tol
        	disp('Minimization complete because optimality tolerance was reached.')
        	break
        end

        k(i) = k(i) + 1;
        k(j) = k(j) - 1;
        b = b + 2*Lp(:, i) - 2*Lp(:,j);
%        b_h = b_h + Hd0L(:,i)-Hd0L(:,j); 

        if opt.Plot
        	E_hist(iter) = (k-k0)'*Lp*(k-k0);
%            rbh = round(b_h)-b_h;
%            rbh = -b_h;
%            E_hist(iter) = E_hist(iter) + rbh'*MM*rbh;
        	m_hist(iter) = m;
        end
    end

    assert(abs(sum(k)-sum(k0)) < 1e-10)

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