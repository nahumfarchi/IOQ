function [alpha_P, beta_P, Lp, out] = ...
    IOQ_highgenus(verts, faces, varargin)
	% Input:
	%	verts - vertices (nv x 3)
	%	faces - faces (nf x 3)
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
    %   'highg_method' <'option1a' | 'option1b' | 'option2'>
    %       The method used for high genus meshes ('round' by default).
    %       'option1a' - optimize E_1 and then use rounding to find beta_P
    %       'option1b' - optimize E_1 and then use cvp to find beta_p (option
    %       1b)
    %       'option2' - set beta_p to be some constant and optimize
    %       both E_1 and E_2 (option 2). Repeat n_alternating times (each
    %       time with the beta_P=round(...)).
    %       'option3' - TODO
    %   'beta_P' <bp|'round'>
    %       The constant value of beta_P used in the alternating method
    %       (zeros by default).
    %   'n_alternating' <n>
	%
	% Output:
	%	

	%%%%% Parse input

	p = inputParser;

	addOptional(p, 'Iterations', 1000);
    addOptional(p, 'NSingularities', []);
    addOptional(p, 'Laplacian', 'conn');
    addOptional(p, 'LaplacianPInv', []);
    addOptional(p, 'Plot', false);
    addOptional(p, 'Tol', 1e-10);
    addOptional(p, 'highg_method', 'option1a');
    addOptional(p, 'beta_P', []);
    addOptional(p, 'n_alternating', 1);
    addOptional(p, 'Verbose', true);
    addOptional(p, 'Histories', false);
    addOptional(p, 'Debug', false);

    parse(p, varargin{:});
    opt = p.Results;
    verb = opt.Verbose;
    mesh = Mesh(verts, faces);
    nv = mesh.nV;
    if nargout > 3 || opt.Plot
        opt.Histories = true;
    end
	E_hist1 = [];
    E_hist2 = [];
    E_round_hist2 = [];
    E_hist = [];
    Emiq_hist = [];
    Emiq_round_hist = [];
    gridsearch_ticks = [];
    m_hist = [];

    %%%%% Setup
    if isempty(opt.LaplacianPInv)
    	if ischar(opt.Laplacian)
    		if verb, disp('Creating L...'); end
    		if strcmpi(opt.Laplacian, 'conn')
    			[d0, ~] = get_exterior_derivatives(mesh);
    			L = d0'*d0;
    		elseif strcmpi(opt.Laplacian, 'cot')
    			L = -cotmatrix(mesh.V, mesh.F);
    		end
    	else
    		L = opt.Laplacian;
    	end

    	if verb, disp('Inverting L...'); end
    	switch opt.InvMethod
        case 'GPUInv'
            Lp = inv(gpuArray(single(L+1/nv))) - 1/nv;
            Lp = double(gather(Lp));
        case 'GPUBlockInv'
            Lp = block_inv_gpu(full(L+1/nv), opt.BlockSize) - 1/nv;
            Lp = double(Lp);
        case 'CholMexInv'
            Lp = invChol_mex(full(L + 1/nv)) - 1/nv;
        otherwise
            error('Unkown InvMethod')
        end
    else
    	Lp = opt.LaplacianPInv;
    end

    if verb, disp('Setting initial singularities...'); end
    alpha_G = get_gaussian_curvature(mesh);
    x0 = (2/pi)*alpha_G;
    if isempty(opt.NSingularities)
    	c = round(abs(sum((x0)))); % = 4 xi
    else
    	c = opt.NSingularities;
    end
    
    alpha_P = zeros(nv, 1);
    n_pos_sing = (c + round(sum(x0))) / 2;
    n_neg_sing = (c - round(sum(x0))) / 2;
    inds_pos = randperm(nv, n_pos_sing);
    inds_neg = randperm(nv, n_neg_sing);
    alpha_P(inds_pos) = 1;
    alpha_P(inds_neg) = -1;
    
    if mesh.genus == 0
        beta_P = [];
    else
        ng2 = 2*mesh.genus;
        [~, K, d0, ~, H] = tcods_gsystem(mesh.V, mesh.F);
        beta_G = K(nv+1:end);
        ne = size(d0, 1);
        %B = H - d0*(d0\H);
        B = (speye(ne) - d0 * Lp * d0') * H;
        %assert(norm(full(d1*d0)) < 1e-10)
        %assert(norm(full(d1*B)) < 1e-10)
        %assert(norm(full(d1*H)) < 1e-10)
        %assert(norm(full(d0'*B)) < 1e-10)
        
        %H = -mesh.H;
        %beta_G = -flatten_generators(mesh);
        
        %if isempty(opt.beta_P)
        %   beta_P = zeros(ng2, 1);
        %elseif ischar(opt.beta_P) && strcmpi(opt.beta_P, 'round')
        %    a = Lp * (alpha_G - (pi/2) * alpha_P);
        %    beta_P = round((beta_G - H'*d0*a) / (pi/2));
        %else
        %    beta_P = opt.beta_P;
        %end
        
        %if ~exist('d0', 'var')
        %    [d0, ~] = get_exterior_derivatives(mesh);
        %end
        
        %y0 = beta_G - (pi/2)*beta_P;        
        %A2 = H'*d0*Lp;
        %M2 = inv(H'*B)*B'*B*inv(H'*B);
    end
    
    %dlp = diag(Lp);
    %b = 2*(Lp*(k - k0));   

    assert(n_pos_sing + n_neg_sing == c)
    assert(abs(sum(alpha_P)-sum(x0)) < 1e-10)

    M1 = Lp;
    %A1 = speye(nv);
    x = alpha_P;
    if opt.Histories || ~strcmp(opt.highg_method, 'option1a')
        %M2 = inv(H'*B)*(B'*B)*inv(H'*B);
        %M2 = B*inv(H'*B);
        %M2 = M2' * M2;
        M2 = B / (H'*B);
        M2 = M2' * M2;
    end
    A2 = H'*d0*Lp;
    %y0 = (2/pi)*(A2*alpha_G - beta_G);
    y0 = A2*(2/pi)*alpha_G - (2/pi)*beta_G;
    if isempty(opt.beta_P)
        y = zeros(ng2, 1);
    elseif ischar(opt.beta_P) && strcmpi(opt.beta_P, 'round')
        y = round(A2*x-y0);
    else
        assert(length(opt.beta_P) == ng2);
        y = opt.beta_P;
    end
    min_val = -inf;
    
    %%%%% Minimize energy
    if verb, disp('Minimizing...'); end
    
    if strcmpi(opt.highg_method, 'option1a')
        assert(ng2 > 0);
        if verb, disp('using option1a'); end
        
        x = gridsearch(M1, x, x0, opt);
        y = round(A2*x - y0);
    elseif strcmpi(opt.highg_method, 'option1b')
        assert(ng2 > 0);
        if verb, disp('cvp'); end
        
        x = gridsearch(M1, x, x0, opt);
        %M2 = B*inv(H'*B);
        target_vec = A2*x-y;
        R = chol(CLLL(M2)'*CLLL(M2));
        target_vec = R*target_vec;
        %lattice_basis = CLLL(chol(M2'*M2));
        %target_vec = lattice_basis \ target_vec;
        
        %R = chol(M2'*M2);
        %target_vec = M2*()
        %lattice_basis = CLLL(chol(M2'*M2));
        %target_vec = lattice_basis*(A2*x-y0);
        
        %B_reduced = CLLL(N2*A2);
        %t = N2*A2*x - N2*y0;
        
        % lattice_solve considers the *rows* to be the basis vectors, so we
        % have to transpose
        %y = lattice_solve(lattice_basis', target_vec); 
        y = babai(R, target_vec);
        y = y(:);
        %fprintf('cvp E2 = %g\n', norm(R*y - target_vec)^2);
        %fprintf('round E2 = %g\n', norm(R*round(A2*x-y0) - target_vec)^2);
        %fprintf('zero E2 = %g\n', norm(R*zeros(ng2, 1) - target_vec)^2);
        %y = lattice_solve(A2*inv(A2'*A2), x-y0);
        %y = round(A2*x - y0);
        %x = lattice_solve(A2', y+y0);
        %y = lattice_solve(A2*inv(A2'*A2), x - inv(A2'*A2)*A2'*y0);
        
        %error('cvp not implemented')
    elseif strcmpi(opt.highg_method, 'option2')
        assert(ng2 > 0);
        if verb, disp('using option2'); end
        
        %M2 = inv(H'*B)*(B'*B)*inv(H'*B);
        %M2 = B*inv(H'*B);
        %M2 = M2' * M2;
        %y = lattice_solve(eye(ng2), A2*x-y0);
        %y = y(:);
        %for ii = 1:opt.n_alternating
            
            x = gridsearch2(M1, x, x0, M2, A2, y, y0, opt);
            %y = lattice_solve(eye(ng2), A2*x-y0);
            %y = y(:);
            %ynew = round(A2*x - y0);
           			
            %if opt.Plot
            %    k = [x; ynew];
            %    res = TCODS(mesh, 'k', k, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
            %end
            
            %if norm(y - ynew) < 1e-10
            %    break
            %end
            
            %y = ynew;
        %end
    elseif strcmpi(opt.highg_method, 'option2_optimized')
        assert(ng2 > 0);
        if verb, disp('using option2_optimized'); end
        
        %M2 = inv(H'*B)*(B'*B)*inv(H'*B);
        %M2 = B*inv(H'*B);
        %M2 = M2' * M2;

        %for ii = 1:opt.n_alternating
            x = gridsearch2_optimized(M1, x, x0, M2, A2, y, y0, opt);
            %y = round(A2*x - y0);
			
            %if opt.Plot
            %    k = [x; y];
            %    res = TCODS(mesh, 'k', k, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
            %end

            %if norm(y - ynew) < 1e-10
            %    break
            %end
            
            %y = ynew;
        %end
    elseif strcmp(opt.highg_method, 'option3')
        assert(ng2 > 0)
        error('Not implemented')
        if verb, disp('using option3'); end
        
        %M2 = inv(H'*B)*(B'*B)*inv(H'*B);

        for ii = 1:opt.n_alternating
            opt.Iterations = 1;
            [alpha_P, E_hist, m_hist] = gridsearch2(M1, x, x0, M2, A2, y0, opt);   
        end
        y = round(A2*x - y0);
    elseif strcmpi(opt.highg_method, 'option3_optimized')
        assert(ng2 > 0);
        if verb, disp('using option3_optimized'); end
        
        %M2 = inv(H'*B)*(B'*B)*inv(H'*B);
        %M2 = B*inv(H'*B);
        %M2 = M2' * M2;
        [x, y, E_hist, m_hist] = gridsearch3_optimized(M1, x, x0, M2, A2, y, y0, opt);
    end
    
    if opt.Plot || opt.Histories
        update_histories(opt.Iterations+1, M1, x, x0, M2, A2, y, y0, min_val)
    end
    
    alpha_P = x;
    beta_P = y;
    
    if opt.Debug
        k = [alpha_P; beta_P];
        res = TCODS(mesh, 'k', k, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
        connection1 = res.connection;
        alpha = alpha_G - (pi/2)*alpha_P;
        beta = beta_G - (pi/2)*beta_P;
        aa = Lp*alpha;
        bb = inv(H'*B)*(beta - H'*d0*aa);
        connection2 = d0*aa + B*bb;
        assert(norm(connection1)-norm(connection2) < 1e-10)
        %assert(norm(connection1-connection2) < 1e-10);

        assert(abs(sum(alpha_P)-sum(x0)) < 1e-10)
    end

    if opt.Plot
        disp('Plotting...')
        xx = 1:length(E_hist1);
        figure(); hold on
        lgd = {};
        if ~isempty(E_hist1) && ~isempty(E_hist2)
            plot(xx, E_hist1+E_hist2)
            lgd{end+1} = 'E1+E2';
        end
        if ~isempty(E_hist1)
            plot(xx, E_hist1)
            lgd{end+1} = 'E1';
        end
        if ~isempty(E_hist2)
            plot(xx, E_hist2)
            lgd{end+1} = 'E2';
        end
        if ~isempty(E_round_hist2)
            plot(xx, E_round_hist2)
            lgd{end+1} = 'E2 w\ y=round(...)';
        end
        if ~isempty(Emiq_hist)
            plot(xx, Emiq_hist)
            lgd{end+1} = 'Emiq';
        end
        if ~isempty(Emiq_round_hist)
            plot(xx, Emiq_round_hist)
            lgd{end+1} = 'Emiq w/ y=round(...)';
        end
        if ~isempty(m_hist)
            plot(xx, m_hist, 'b')
            lgd{end+1} = 'm';
        end
        hold off
        legend(lgd{:})
        xlabel('Iteration')
        title('IOQ (gpu)')
    end
    
    if opt.Histories
        out.E_hist1 = E_hist1;
        out.E_hist2 = E_hist2;
        out.E_round_hist2 = E_round_hist2;
        out.E_hist = E_hist;
        out.Emiq_hist = Emiq_hist;
        out.Emiq_round_hist = Emiq_round_hist;
        out.gridsearch_ticks = gridsearch_ticks;
        out.m_hist = m_hist;
    end

    if verb, fprintf('Done.\n\n'); end

    
    
    
    
    
    
    
    
    
    
function x = gridsearch(M, x, x0, opt)
    % Minmize (Ax - x0)^T M (Ax - x0)

    M_tilde = M;
    b = 2*M*(x - x0);
    dm = diag(M_tilde);
    for iter = 1:opt.Iterations       
        if opt.Histories
            update_histories(iter, M1, x, x0, M2, A2, y, y0, min_val);
        end

        if verb && mod(iter, 10)==0
            fprintf('iter, m : %d, %g\n', iter, min_val);
        end
        
        [min_val, i] = min(bsxfun(@plus, dm+b, (dm-b)') - 2*M_tilde);
        [min_val, j] = min(min_val);
        i = i(j);
        
        
        
        if abs(min_val(1)) < opt.Tol
            disp('Minimization complete because optimality tolerance was reached.')
            break
        end
        
        x(i) = x(i) + 1;
        x(j) = x(j) - 1;
        b = b + 2*M_tilde(:,i) - 2*M_tilde(:,j);
    end
    
    if iter == opt.Iterations
        fprintf('Optimization terminated because the iterations budget %d was reached\n', opt.Iterations);
        fprintf('min_val = %g\n', min_val);
    end

end

    
function x = gridsearch2(M1, x, x0, M2, A2, y, y0, opt)
    % Minimize (A1*x - x0)^T M1 (A1*x - x0) + (A2*y - y0)^T M2 (A2*y - y0)

    y0_tilde = y + y0;
    M1_tilde = M1;
    b1 = 2*M1*(x - x0);
    dm1 = diag(M1_tilde);
    M2_tilde = A2'*M2*A2;
    b2 = 2*A2'*M2'*(A2*x - y0_tilde);
    dm2 = diag(M2_tilde);
    for iter = 1:opt.Iterations
        if opt.Histories
            update_histories(iter, M1, x, x0, M2, A2, y, y0, min_val);
        end
        
        [min_val, i] = min(bsxfun(@plus, dm1+b1, (dm1-b1)') - 2*M1_tilde + ...
                           (bsxfun(@plus, dm2+b2, (dm2-b2)') - 2*M2_tilde));
        [min_val, j] = min(min_val);
        i = i(j);
        
        if verb && mod(iter, 10)==0
            fprintf('iter, m : %d, %g\n', iter, min_val);
        end
        
        if abs(min_val(1)) < opt.Tol
            disp('Minimization complete because optimality tolerance was reached.')
            break
        end
        
        x(i) = x(i) + 1;
        x(j) = x(j) - 1;
        b1 = b1 + 2*M1_tilde(:,i) - 2*M1_tilde(:,j);
        b2 = b2 + 2*M2_tilde(:,i) - 2*M2_tilde(:,j);
    end
    
    if iter == opt.Iterations
        disp('Minimization complete because the max number of iterations was reached.')
    end
    
end

function x = gridsearch2_optimized(M1, x, x0, M2, A2, y, y0, opt)
    % Minimize (x-x0)^T M1 (x-x0) + (A2 x - y0_tilde)^T M2 (A2 x -
    % y0_tilde)
    
    y0_tilde = y + y0;
    F = M1 + A2'*M2*A2;
    df = diag(F);
    f = M1*x0 + A2'*M2*y0_tilde;
    g = 2*F'*x - 2*f;
    for iter = 1:opt.Iterations
        [min_val, i] = min(bsxfun(@plus, df+g, (df-g)') - 2*F);
        [min_val, j] = min(min_val);
        i = i(j);
        
        if verb && mod(iter, 10)==0
            fprintf('iter, m : %d, %g\n', iter, min_val);
        end
        
        if abs(min_val(1)) < opt.Tol
            disp('Minimization complete because optimality tolerance was reached.')
            break
        end
        
        x(i) = x(i) + 1;
        x(j) = x(j) - 1;
        g = g + 2*F(:, i) - 2*F(:, j);
        
        if opt.Histories
            update_histories(iter, M1, x, x0, M2, A2, y, y0, min_val)
        end
    end
    
    if iter == opt.Iterations
        disp('Minimization complete because the max number of iterations was reached.')
    end
end
    
function [x, y, E_hist, m_hist] = gridsearch3_optimized(M1, x, x0, M2, A2, y, y0, opt)
    % Minimize (x-x0)^T M1 (x-x0) + (A2 x - y0_tilde)^T M2 (A2 x -
    % y0_tilde)
    
    y0_tilde = y + y0;
    F = M1 + A2'*M2*A2;
    df = diag(F);
    f = M1*x0 + A2'*M2*y0_tilde;
    g = 2*F'*x - 2*f;
    for iter = 1:opt.Iterations
        [min_val, i] = min(bsxfun(@plus, df+g, (df-g)') - 2*F);
        [min_val, j] = min(min_val);
        i = i(j);
        
        if verb && mod(iter, 10)==0
            fprintf('iter, m : %d, %g\n', iter, min_val);
        end
        
        if abs(min_val(1)) < opt.Tol
            disp('Minimization complete because optimality tolerance was reached.')
            break
        end
        
        x(i) = x(i) + 1;
        x(j) = x(j) - 1;
        y = round(A2*x - y0);
        y0_tilde = y + y0;
        f = M1*x0 + A2'*M2*y0_tilde;
        g = 2*F'*x - 2*f;
        %g = g + 2*F(:, i) - 2*F(:, j);
        
        if opt.Histories
            update_histories(iter, M1, x, x0, M2, A2, y, y0, min_val)
        end
    end
    
    if iter == opt.Iterations
        disp('Minimization complete because the max number of iterations was reached.')
    end 
end

function update_histories(iter, M1, x, x0, M2, A2, y, y0, min_val)
    y0_tilde = y + y0;
    E_hist1(end+1) = gather((x-x0)'*M1*(x-x0));
    E_hist2(end+1) = gather((A2*x-y0_tilde)'*M2*(A2*x-y0_tilde));
    E_hist(end+1) = E_hist1(end) + E_hist2(end);
    k = [x; y];
    res = TCODS(mesh, 'k', k, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
    Emiq_hist(end+1) = res.miq_energy;
    yround = round(A2*x-y0);
    k = [x; yround];
    res = TCODS(mesh, 'k', k, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
    Emiq_round_hist(end+1) = res.miq_energy;
    E_round_hist2(end+1) = gather((A2*x-yround-y0)'*M2*(A2*x-yround-y0));
    m_hist(end+1) = gather(min_val);
    if iter == 1
        gridsearch_ticks(end+1) = 1;
    else
        gridsearch_ticks(end+1) = 0;
    end
    
    if opt.Debug
        alpha_ = (2/pi)*alpha_G - x;
        beta_ = (2/pi)*beta_G - y;
        aa_ = Lp*alpha_;
        bb_ = inv(H'*B)*(beta_ - H'*d0*aa_);
        assert(norm(d0*aa_)^2 - E_hist1(end) < 1e-10)
        assert(norm(B*bb_)^2 - E_hist2(end) < 1e-10)
        
        tmp = (H'*d0*Lp*x - (y - (2/pi)*beta_G + H'*d0*Lp*(2/pi)*alpha_G));
        assert(tmp'*M2*tmp - E_hist2(end) < 1e-10)
        tmp = (2/pi)*alpha_G - x;
        assert(tmp'*Lp*tmp - E_hist1(end) < 1e-10)
        
        assert(Emiq_hist(end) - E_hist(end)*(pi/2)^2 < 1e-10)
        
        % |x_tcods|^2 = (pi/2)^2 * |d0 a + Bb|^2
        assert(Emiq_hist(end) - (pi/2)^2*norm(d0*aa_ + B*bb_)^2 < 1e-10)
        % |x_tcods|^2 = (pi/2)^2 * ( |d0 a|^2 + |Bb|^2 )
        assert(Emiq_hist(end) - (pi/2)^2*(norm(d0*aa_)^2 + norm(B*bb_)^2) < 1e-10)
    end
end

end
  
% function [x, E_hist, m_hist] = gridsearch3(M1, x, x0, M2, A2, y, y0, opt)
%     % Minimize (x-x0)^T M1 (x-x0) + (A2 x - y0_tilde)^T M2 (A2 x -
%     % y0_tilde)
%     if opt.Plot
%         E_hist = zeros(opt.Iterations, 1);
%         m_hist = zeros(opt.Iterations, 1);
%     else
%         E_hist = [];
%         m_hist = [];
%     end
%     
%     y0_tilde = y + y0;
%     F = A1'*M1*A1 + A2'*M2*A2;
%     df = diag(F);
%     f = M1*x0 + A2*M2*y0_tilde;
%     for iter = 1:opt.Iterations
%         [min_val, i] = min(bsxfun(@plus, g-df, (g-df)'), - 2*F);
%         [min_val, j] = min(min_val);
%         i = i(j);
%         
%         if abs(min_val(1)) < opt.Tol
%             disp('Minimization complete because optimality tolerance was reached.')
%             break
%         end
%         
%         x(i) = x(i) + 1;
%         x(j) = x(j) - 1;
%         g = g + 2*F(:, i) - 2*F(:, j);
%         
%         if opt.Plot
%             E_hist(iter) = (A1*x-x0)'*M1*(A1*x-x0) + (A2*x-y0_tilde)'*M2*(A2*x-y0_tilde);
%             m_hist(iter) = min_val;
%         end
%     end
%     
%     if iter == opt.Iterations
%         disp('Minimization complete because the max number of iterations was reached.')
%     end
%     
%     if opt.Plot
%         E_hist = E_hist(1:iter-1);
%         m_hist = m_hist(1:iter-1);
%     end
    
% function [x, E_hist, m_hist] = gridsearch3(M1, A1, x, x0, M2, A2, y0, opt)
%     % Minimize (A1*x - x0)^T M1 (A1*x - x0) + (A2*y - y0)^T M2 (A2*y - y0)
%     if opt.Plot
%         E_hist = zeros(opt.Iterations, 1);
%         m_hist = zeros(opt.Iterations, 1);
%     else
%         E_hist = [];
%         m_hist = [];
%     end
%     
%     %y0 = (2/pi)*(H'*d0*Lp*alpha_G+beta_G);
%     y = round(A2*x - y0);
%     y0_tilde = y0 - y;
%     
%     M1_tilde = A1'*M1*A1;
%     b1 = 2*A1'*M1*(A1*x - x0);
%     dm1 = diag(M1_tilde);
%     M2_tilde = A2'*M2*A2;
%     b2 = 2*A2'*M2'*(A2*x - y0_tilde);
%     dm2 = diag(M2_tilde);
%     for iter = 1:opt.Iterations
%         [min_val, i] = min(bsxfun(@plus, dm1+b1, (dm1-b1)') - 2*M1_tilde + ...
%                            bsxfun(@plus, dm2+b2, (dm2-b2)') - 2*M2_tilde);
%         [min_val, j] = min(min_val);
%         i = i(j);
%         
%         if abs(min_val(1)) < opt.Tol
%             disp('Minimization complete because optimality tolerance was reached.')
%             break
%         end
%         
%         x(i) = x(i) + 1;
%         x(j) = x(j) - 1;
%         b1 = b1 + 2*M1_tilde(:,i) - 2*M1_tilde(:,j);
%         b2 = b2 + 2*M2_tilde(:,i) - 2*M2_tilde(:,j);
%         
%         y = round(A2*x - y0);
%         y0_tilde = y0 - y;
%         
%         if opt.Plot
%             E_hist(iter) = (A1*x-x0)'*M1*(A1*x-x0) + (A2*x-y0_tilde)'*M2*(A2*x-y0_tilde);
%             m_hist(iter) = min_val;
%         end
%     end
%     
%     if iter == opt.Iterations
%         disp('Minimization complete because the max number of iterations was reached.')
%     end
%     
%     if opt.Plot
%         E_hist = E_hist(1:iter-1);
%         m_hist = m_hist(1:iter-1);
%     end



