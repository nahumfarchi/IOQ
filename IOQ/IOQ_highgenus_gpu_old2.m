function [alpha_P, beta_P, elapsed_total, Lp, out] = ...
    IOQ_highgenus_gpu(verts, faces, varargin)
	% Input:
	%	verts - (nv x 3)
	%	faces - (nf x 3)
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
    %       'option2_optimized'
    %       'option3' - TODO
    %   'beta_P' <bp|'round'>
    %       The constant value of beta_P used in the alternating method
    %       (zeros by default).
    %   'n_alternating' <n>
	%
	% Output:
	%   TODO	
    %
    % Example:
    %   m = Mesh(path_to_mesh);
    %   [alpha_p, beta_p] = IOQ_highgenus_gpu(m.V, m.F, ...)
    %   k = [alpha_p; beta_p];
    %   res = TCODS(m, 'k', k, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);

	%%%%% Parse input

	parser = inputParser;

	addOptional(parser, 'Iterations', 1000);
    addOptional(parser, 'NSingularities', []);
    addOptional(parser, 'Laplacian', 'conn');
    addOptional(parser, 'LaplacianPInv', []);
    addOptional(parser, 'Plot', false);
    addOptional(parser, 'Tol', 1e-6);
    addOptional(parser, 'highg_method', 'option1a');
    addOptional(parser, 'beta_P', []);
    addOptional(parser, 'n_alternating', 1);
    addOptional(parser, 'BlockSize', 20000);
    addOptional(parser, 'UseGPU', true);
    addOptional(parser, 'bsx', false);
    addOptional(parser, 'Kernel', true);
    addOptional(parser, 'Verbose', true);
    addOptional(parser, 'Histories', false);
    addOptional(parser, 'Debug', false);
    addOptional(parser, 'InvMethod', 'GPUInv');
    addOptional(parser, 'JLEps', 0.5);
    addOptional(parser, 'JLFac', 24);
    addOptional(parser, 'SampleResistance', false);
    addOptional(parser, 'NSamples', 1000);
    addOptional(parser, 'Mesh', []);
    addOptional(parser, 'KernelEntry', []);
    addOptional(parser, 'Colamd', false);
    addOptional(parser, 'CholGPU', false);
    addOptional(parser, 'ZtildeOneBS', false);
    addOptional(parser, 'block_gs', false);

    parse(parser, varargin{:});
    opt = parser.Results;
    opt.Kernel = ~opt.bsx;
    verb = opt.Verbose;
    if isempty(opt.Mesh)
        mesh = Mesh(verts, faces);
    else
        mesh = opt.Mesh;
    end
    nv = mesh.nV; ne = mesh.nE;
    if nargout > 4 || opt.Plot
        opt.Histories = true;
    end
    E_hist1 = []; E_hist2 = []; E_round_hist2 = []; E_hist = [];
    Emiq_hist = []; Emiq_round_hist = []; gridsearch_ticks = [];
    m_hist = [];
    elapsed_total = 0;
    if isempty(opt.highg_method)
        opt.highg_method = 'option1a';
    end
    switch opt.InvMethod
        case 'ApproxResistance'
            opt.LaplacianPInv = []; % make sure this is ignored
            if isempty(opt.KernelEntry)
                opt.KernelEntry = 'reduce_cols_R_symmetric';
            end
        case 'NoInv'
            opt.LaplacianPInv = [];
            opt.KernelEntry = 'reduce_cols'; % TODO ?
            Lp = [];
        otherwise
            opt.KernelEntry = 'reduce_cols';
    end

    %%%%% Setup
    
    % Setup the gpu kernel
    tic
    if opt.UseGPU && (opt.Kernel || opt.bsx)
        kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', opt.KernelEntry);
        num_elem = nv^2;
        kernel.ThreadBlockSize = [1, kernel.MaxThreadsPerBlock, 1];
        kernel.GridSize = [1, ceil(nv / kernel.MaxThreadsPerBlock)];
        out_min = single(zeros(1, nv, 'gpuArray'));
        out_r = uint32(zeros(1, nv, 'gpuArray'));
        gd = gpuDevice();
    else
        opt.invMethod = 'CholMexInv';
        opt.Kernel = false;
        gd = [];
    end
    
    if verb, disp('Setting initial singularities...'); end
    
    alpha_G = get_gaussian_curvature(mesh);
    genus = round(1 - sum(alpha_G)/(4*pi));
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
    
    if ~isempty(gd), wait(gd); end; elapsed_total = elapsed_total + toc;
    
    if isempty(opt.LaplacianPInv)
        tic
    	if ischar(opt.Laplacian)
    		if verb, disp('Creating L...'); end
    		if strcmpi(opt.Laplacian, 'conn')
    			d0 = get_exterior_derivatives(mesh);
                if strcmpi(opt.InvMethod, 'ApproxResistance') && opt.Colamd
                    p = colamd(d0);
                    pt(p) = 1:length(p);
                    d0p = d0(:,p);
                    Lperm = d0p' * d0p;
                end
    			L = d0'*d0;
    		elseif strcmpi(opt.Laplacian, 'cot')
    			L = -cotmatrix(mesh.V, mesh.F);
    		end
        else
    		L = opt.Laplacian;
    	end

    	if verb, disp('Inverting L...'); end
    	%Lp = invChol_mex(full(L + 1/nv)) - 1/nv;
        
        switch opt.InvMethod
        case 'GPUInv'
            Lp = inv(gpuArray(single(L+1/nv))) - 1/nv;
            Lp = double(gather(Lp));
        case 'GPUBlockInv'
            Lp = block_inv_gpu(full(L+1/nv), opt.BlockSize, true) - 1/nv;
            Lp = double(Lp);
        case 'ApproxResistance'
            F = factorize(L);
            [Ztilde, Rtilde] = resistance_distance();
        case 'CholMexInv'
            Lp = invChol_mex(full(L + 1/nv)) - 1/nv;
        case 'NoInv'
            %Lchol = ichol(Lperm,struct('type','ict','droptol',1e-4,'michol','off'));
            [Lchol, ~, S] = chol(L);
        otherwise
            error('Unkown InvMethod')
        end
        if ~isempty(gd), wait(gd); end; elapsed_total = elapsed_total + toc;
    else
    	Lp = opt.LaplacianPInv;
    end
    
    tic
    if genus == 0
        beta_P = [];
    else
        ng2 = 2*genus;
        if ~exist('d0', 'var')
            d0 = get_exterior_derivatives(mesh);
        end
        %Kg = get_gaussian_curvature(mesh);
        %z = -flatten_generators(mesh);
        H = -mesh.H;
        z = wrapToPi(generator_angle_defects(mesh));
        %z = generator_angle_defects(mesh);
        K = [alpha_G; z];
        beta_G = K(nv+1:end);
        %ne = size(d0, 1);
        %B = (speye(ne) - d0 * Lp * d0') * H;
        %if exist('Lp', 'var')
        %    B = H - d0 * (Lp * (d0' * H)); % much faster this way
        %else
        %    B = H - d0 *
        %    %B = H - d0*(d0\H);
        %end
        if ~exist('F', 'var')
            F = factorize(L);
        end
        B = H - d0 * (F \ (d0' * H));
         
        %beta_G = -flatten_generators(mesh);
        %beta_G = generator_angle_defects(mesh);
        
        %if isempty(opt.beta_P)
        %   beta_P = zeros(ng2, 1);
        %elseif ischar(opt.beta_P) && strcmpi(opt.beta_P, 'round')
        %    a = Lp * (alpha_G - (pi/2) * alpha_P);
        %    beta_P = round((beta_G - H'*d0*a) / (pi/2));
        %else
        %    beta_P = opt.beta_P;
        %end
        
        %y0 = beta_G - (pi/2)*beta_P;        
        %A2 = H'*d0*Lp;
        %M2 = inv(H'*B)*B'*B*inv(H'*B);
    end
    if ~isempty(gd), wait(gd); end; elapsed_total = elapsed_total + toc;
    
    %dlp = diag(Lp);
    %b = 2*(Lp*(k - k0));   

    assert(n_pos_sing + n_neg_sing == c)
    assert(abs(sum(alpha_P)-sum(x0)) < 1e-10)
    
    
    
    %Lp = gpuArray(single(Lp));
    %dlp = gpuArray(single(dlp));
    %b = gpuArray(single(b));
    %k = gpuArray(single(k));
    %k0 = gpuArray(single(k0));

    %M1 = Lp;
    %A1 = speye(nv);
    
    
    %%%%% Minimize energy
    if verb, disp('Minimizing...'); end
    
    min_val = inf;
    x = alpha_P;
    if genus == 0
        if ~strcmpi(opt.highg_method, 'genus0')
            error(['highg_method is ', opt.highg_method, ' but genus is 0'])
        end
        if opt.Histories && (~exist('Lp', 'var') || size(Lp,1)~=nv)
            if opt.UseGPU
                Lp = inv(gpuArray(single(L+1/nv))) - 1/nv;
            else
                Lp = invChol_mex(full(L+1/nv)) - 1/nv;
            end
        end
        tic
        if opt.block_gs
            alpha_P = gridsearch_blocks(Lp, x, x0, opt);
            alpha_P = gather(alpha_P);
        else
            switch opt.InvMethod
                case 'ApproxResistance'
                    alpha_P = gridsearch_approx_res(Rtilde, F, Ztilde, x, x0, opt);
                    alpha_P = gather(alpha_P);
                    %if opt.Colamd
                    %    alpha_P = alpha_P(pt);
                    %end
                case 'NoInv'
                    alpha_P = gridsearch_noinv(L, Lchol, S, x, x0, opt);
                otherwise
                    alpha_P = gridsearch(Lp, x, x0, opt);
                    alpha_P = gather(alpha_P);
            end
        end
        beta_P = [];
        if ~isempty(gd), wait(gd); end; elapsed_total = elapsed_total + toc;
    else
        tic
        if opt.Histories || (~strcmp(opt.highg_method, 'option1a'))
            M2 = B / (H'*B);
            M2 = M2' * M2;
        end
        A2 = H'*d0*Lp;
        y0 = (2/pi)*(A2*alpha_G - beta_G);
        if isempty(opt.beta_P)
            y = zeros(ng2, 1);
        elseif ischar(opt.beta_P) && strcmpi(opt.beta_P, 'round')
            y = round(A2*x-y0);
        else
            assert(length(opt.beta_P) == ng2);
            y = opt.beta_P;
        end
        
        switch opt.highg_method
            case 'option1a'       
        %if strcmpi(opt.highg_method, 'option1a')
                assert(ng2 > 0);
                if verb, disp('using option1a'); end
            
                x = gridsearch(Lp, x, x0, opt);
                y = round(A2*x - y0);
            case 'option1b'
        %elseif strcmpi(opt.highg_method, 'option1b')
                assert(ng2 > 0);
                if verb, disp('cvp'); end

                x = gridsearch(Lp, x, x0, opt);
                %M2 = B*inv(H'*B);
                target_vec = A2*x-y;
                R = chol(CLLL(M2)'*CLLL(M2));
                target_vec = R*target_vec;
                y = babai(R, target_vec);
                y = y(:);
                %y = lattice_solve(eye(ng2), A2*x-y0);
                %y = y(:);
                %y = lattice_solve(A2*inv(A2'*A2), x-y0);
                %y = round(A2*x - y0);
                %x = lattice_solve(A2', y+y0);
                %y = lattice_solve(A2*inv(A2'*A2), x - inv(A2'*A2)*A2'*y0);

                %error('cvp not implemented')
            case 'option2'
        %elseif strcmpi(opt.highg_method, 'option2')
                assert(ng2 > 0);
                if verb, disp('using option2'); end

                %M2 = inv(H'*B)*(B'*B)*inv(H'*B);
                %M2 = B*inv(H'*B);
                %M2 = M2' * M2;
                %y = lattice_solve(eye(ng2), A2*x-y0);
                %y = y(:);
                %for ii = 1:opt.n_alternating

                    x = gridsearch2(Lp, x, x0, M2, A2, y, y0, opt);
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
            case 'option2_optimized'
        %elseif strcmpi(opt.highg_method, 'option2_optimized')
                assert(ng2 > 0);
                if verb, disp('using option2_optimized'); end

                %M2 = inv(H'*B)*(B'*B)*inv(H'*B);
                %M2 = B*inv(H'*B);
                %M2 = M2' * M2;

                %for ii = 1:opt.n_alternating
                    x = gridsearch2_optimized(Lp, x, x0, M2, A2, y, y0, opt);
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
            case 'option3'
        %elseif strcmp(opt.highg_method, 'option3')
                assert(ng2 > 0)
                error('Not implemented')
                if verb, disp('using option3'); end

                M2 = inv(H'*B)*(B'*B)*inv(H'*B);

                for ii = 1:opt.n_alternating
                    opt.Iterations = 1;
                    alpha_P = gridsearch2(Lp, x, x0, M2, A2, y0, opt);   
                end
                y = round(A*x - y0);
            case 'option3_optimized'
        %elseif strcmpi(opt.highg_method, 'option3_optimized')
                assert(ng2 > 0);
                if verb, disp('using option3_optimized'); end

                %M2 = inv(H'*B)*(B'*B)*inv(H'*B);
                %M2 = B*inv(H'*B);
                %M2 = M2' * M2;
                [x, y] = gridsearch3_optimized(Lp, x, x0, M2, A2, y, y0, opt);
            case 'genus0'
                error('genus0 method was given but mesh is high genus')
            otherwise
                error('Unkown highg_method was given')
        end
        
        alpha_P = gather(double(x));
        beta_P = gather(double(y));
        
        if ~isempty(gd), wait(gd); end; elapsed_total = elapsed_total + toc;

        if opt.Plot || opt.Histories
            update_histories(opt.Iterations+1, Lp, x, x0, M2, A2, y, y0, min_val)
        end

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
            if ~isempty(Emiq_hist)
                plot(xx, Emiq_hist)
                lgd{end+1} = 'Emiq';
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
    end

    if opt.Histories
        out.E_hist1 = E_hist1;
        out.E_hist2 = E_hist2;
        out.E_hist = E_hist;
        out.Emiq_hist = Emiq_hist;
        out.gridsearch_ticks = gridsearch_ticks;
        out.m_hist = m_hist;
    end

    if verb, fprintf('Done.\n\n'); end
    
function x = gridsearch_approx_res(Rtilde, F, Ztilde, x, x0, opt)
    if opt.UseGPU && ~opt.SampleResistance
       Rtilde = single(gpuArray(Rtilde));
    end
    %Fperm = factorize(Lperm);

    tol = 1e-6;
    maxit = 100;

    i_hist = [];
    j_hist = [];
    btilde_cpu = [];
    
    for iter = 1:opt.Iterations
        if opt.Histories
            E_hist1(end+1) = gather((x-x0)'*Lp*(x-x0));
        end

        if opt.Verbose && mod(iter, 10) == 0
            fprintf('iter, m : %d, %g\n', iter, min_val);
        end

        % Not sure why but doing this on the GPU is slower
        rhs = 2*(x-x0);
        if opt.Colamd
            rhs = rhs(p);
        end
        %btilde_cpu = Lchol' \ (Lchol \ rhs);
        %[btilde_cpu, flag] = pcg(Lperm, rhs, 1e-6, 100, Lchol, Lchol', gather(btilde_cpu));
        btilde_cpu = F \ rhs;
        if opt.Colamd
            btilde_cpu = btilde_cpu(pt);
        end
        %if opt.Colamd
        %    btilde = pcg(Lperm, rhs, tol, maxit, Lchol, Lchol', double(gather(btilde)));
        %    btilde = btilde(pt);
        %else
        %    btilde = pcg(L, rhs, tol, maxit, Lchol, Lchol', double(gather(btilde)));
        %end
        if opt.UseGPU
            btilde = single(gpuArray(btilde_cpu));
        end
        if opt.SampleResistance
            rows = randi(nv, opt.NSamples, 1);
            cols = randi(nv, opt.NSamples, 1);
            min_val = inf;
            for samp = 1:opt.NSamples
                r = rows(samp);
                c = cols(samp);
                if Rtilde(r, c) == 0
                    Rtilde(r, c) = gather(norm(Ztilde(:, r) - Ztilde(:, c))^2);
                end
                m = btilde(r) - btilde(c) + Rtilde(r, c);
                if m < min_val
                    min_val = m;
                    i = r;
                    j = c;
                end
            end
        elseif opt.Kernel && opt.UseGPU
            [out_min, out_r] = feval(kernel, out_min, out_r, Rtilde, btilde, nv, nv, num_elem);
            %out_r = out_r;
            [min_val, j] = min(out_min);
            i = out_r(j) + 1;
        else
            %[min_val, i] = min(bsxfun(@plus, btilde, -btilde') + Rtilde);
            %[min_val, i] = min(gather(bsxfun(@plus, btilde, -btilde')) + Rtilde)
            [min_val, i] = min(bsxfun(@plus, btilde, -btilde') + squareform(Rtilde));
            [min_val, j] = min(min_val);
            i = i(j);
            fprintf('%d, %d\n', i, j);
        end

        i_hist(iter) = gather(i);
        j_hist(iter) = gather(j);
        if length(i_hist) > 2 && length(j_hist) > 2 && i_hist(iter-2) == i && j_hist(iter-2) == j
            disp('Cycle detected. Stopping minimization');
            disp(['m : ', num2str(min_val)])
            break;
        end
        
        if abs(min_val(1)) < opt.Tol
            disp('Minimization complete because optimality tolerance was reached.')
            disp(['m : ', num2str(min_val)])
            break
        end
    
        if i == j
            disp('Minimization complete because i==j');
            break
        end
    
        x(i) = x(i) + 1;
        x(j) = x(j) - 1;
    end

    if iter == opt.Iterations
        fprintf('Optimization terminated because the iterations budget %d was reached\n', opt.Iterations);
        fprintf('min_val = %g\n', min_val);
    end
end

function x = gridsearch_blocks(M, x, x0, opt)
    sz = 10000;
    M_cell = mat2tiles(M, nv, sz);
    if ~isa(M, 'single')
        M = single(M);
    end
    b = gpuArray(2*M*(x - x0));
    dm = gpuArray(single(diag(M)));
    b_cell = mat2tiles(b, sz, 1);
    dm_cell = mat2tiles(dm, sz, 1);
    %M_blk = zeros(nv, sz, 'like', b);
    
    for iter = 1:opt.Iterations
        min_val = inf;
        best_i = nan;
        best_j = nan;
        
        for blk = 1:numel(M_cell)
            
            %cols = (blk-1)*sz + 1 : min(blk*sz, size(;
            %M_blk = gpuArray(M_cell{blk});
            M_cell{blk} = gpuArray(M_cell{blk});
            b_blk = b_cell{blk};
            dm_blk = dm_cell{blk};
            [min_val_blk, i] = min(bsxfun(@plus, dm+b, (dm_blk-b_blk)') - 2*M_cell{blk});
            [min_val_blk, j] = min(min_val_blk);
            i = i(j);
            M_cell{blk} = gather(M_cell{blk});
            
%             [out_min, out_r] = feval(kernel, out_min, out_r, M, dm, b, nv, nv, num_elem);
%             out_r = out_r + 1;
%             [min_val, j] = min(out_min);
%             i = out_r(j);
            
            if min_val_blk < min_val
                min_val = min_val_blk;
                best_i = i;
                best_j = j + (blk-1)*sz;
            end
        end
        
        if abs(min_val) < opt.Tol
            disp('Minimization complete because optimality tolerance was reached.')
            disp(['m : ', num2str(min_val)])
            break
        end
        
        x(best_i) = x(best_i) + 1;
        x(best_j) = x(best_j) - 1;
        b = b + 2*M(:,best_i) - 2*M(:,best_j);
        b_cell = mat2tiles(b, sz, 1);
    end
end

function x = gridsearch_noinv(L, Lchol, S, x, x0, opt)
    % M - Lp
    % x - alpha_p
    % x0 - 
    
    %Lp = invChol_mex(full(L + 1/nv)) - 1/nv;
    
    min_val = inf;
    ener = inf;
    best_i = nan;
    best_j = nan;
    b = x - x0;
    %b = b(p);
    
    % Create a matrix with all possible eij's as columns (each eij is a
    % vector with e(i) = 1 and e(j) = -1).
    
    %J = combnk(1:nv, 2);
    %cols = size(J, 1);
    %I = [1:cols, 1:cols]';
    %J = J(:);
    %V = [ones(cols, 1); -ones(cols, 1)];
    %B0 = sparse(I, J, V, cols, nv, 2*cols)';
    %B0 = [B0, -B0];
    J = combinator(nv, 2, 'p');
    cols = size(J, 1);
    I = [1:cols, 1:cols]';
    V = [ones(cols, 1); -ones(cols, 1)];
    B0 = sparse(I, J, V, cols, nv, 2*cols)';
    F = factorize(L);
    for iter = 1:opt.Iterations
        clear sol;
        %RHS = repmat(b, 1, 2*cols) + B0; % TODO improve memory consumption
        
        %fprintf('%d, %.4g, %.4g\n', iter, ener, b'*Lp*b)
        
        %RHS_cell = mat2tiles(RHS, nv, 128);
        %parfor blk = 1:numel(RHS_cell)
        %    sol{blk} = S * (Lchol \ (Lchol' \ (S' * RHS_cell{blk})));
        %end
        %X = [sol{:}];
        
        %for i = 1:nv
        %    disp(i)
        %    for j = 1:nv
        %        sol = S * (Lchol \ (Lchol' \ (S' * b)));
        %    end
        %end
        
        %X = S * (Lchol \ (Lchol' \ (S' * RHS)));
        X = F \ B0;
        
        fprintf('err : %.4g\n', norm(L*X-RHS, 'fro') / norm(RHS, 'fro'));
        
        %energies = diag(RHS' * X);
        energies = sum(RHS' .* X', 2);
        [ener, col] = min(energies);
        %if abs(ener - min_val) > opt.Tol
        if ener - min_val < -opt.Tol
            min_val = ener;
            eij = B0(:, col);
            x = x + eij;
            b = b + eij;
        else
            disp('Minimization complete because optimality tolerance was reached.')
            break
        end
    end
    
    % ichol + amd + pcg
%     for iter = 1:opt.Iterations
%         [z0, flag] = pcg(L, b, 1e-6, 100, Lchol, Lchol');
%         
%         for i = 1:nv            
%             b(i) = b(i) + 1;
%             [zi, flag] = pcg(L, b, 1e-6, 100, Lchol, Lchol', z0);
%             
%             for j = 1:nv
%                 b(j) = b(j) - 1;
%                 [zij, flag] = pcg(L, b, 1e-6, 100, Lchol, Lchol', zi);
%                 energy = b' * zij;
%                 if energy < min_val
%                     min_val = energy;
%                     best_i = i;
%                     best_j = j;
%                 end
%                 b(j) = b(j) + 1;
%             end
%             
%             b(i) = b(i) - 1;
%         end
%         
%         b(best_i) = b(best_i) + 1;
%         b(best_j) = b(best_j) - 1;
%         min_val = inf;
%     end
%     
%     b = b(pt);
%     x = b;
end

function x = gridsearch(M, x, x0, opt)
    % Minmize (Ax - x0)^T M (Ax - x0)
    % M has to be symmetric (this is not checked)

    if opt.UseGPU
        try
            % Calculate b and dm on the cpu first (uses less gpu memory)
            b = gpuArray(single(2*M*(x - x0)));
            dm = gpuArray(single(diag(M)));
            M = gpuArray(single(M));
        catch
            % Probably out of memory, so try to pass only the lower
            % triangular part of M (it's symmetric).
            %M = gather(M);
            %b = gpuArray(single(2*M*(x - x0)));
            %dm = gpuArray(single(diag(M)));
            % we are now passing the resistance to the kernel
            dm = gather(dm);
            M = bsxfun(@plus, dm, dm) - 2*M;  
            M = tril(M, -1);
            % Just the lower part (not including the diagonal)
            M = M(M > 0);
            M = gpuArray(single(M));
            dm = gpuArray(dm);
                
            % Setup the symmetric kernel
            opt.KernelEntry = 'reduce_cols_R_symmetric';
            kernel = parallel.gpu.CUDAKernel('gridsearch_inner.ptx', 'gridsearch_inner.cu', opt.KernelEntry);
            num_elem = nv^2;
            kernel.ThreadBlockSize = [1, kernel.MaxThreadsPerBlock, 1];
            kernel.GridSize = [1, ceil(nv / kernel.MaxThreadsPerBlock)];
        end
    else
        b = 2*M*(x - x0);
        dm = diag(M);
    end
    
    for iter = 1:opt.Iterations
        if opt.Histories
            if genus > 0
                update_histories(iter, M, x, x0, M2, A2, y, y0, min_val);
            else
                E_hist1(end+1) = gather((x-x0)'*M*(x-x0));
            end
        end

        if opt.Verbose && mod(iter, 10) == 0
            fprintf('iter, m : %d, %g\n', iter, min_val);
        end

        if (opt.UseGPU && opt.bsx) || ~opt.UseGPU
            [min_val, i] = min(bsxfun(@plus, dm+b, (dm-b)') - 2*M);
            [min_val, j] = min(min_val);
            i = i(j);
        elseif opt.Kernel
            switch opt.KernelEntry
                case 'reduce_cols'
                    [out_min, out_r] = feval(kernel, out_min, out_r, M, dm, b, nv, nv, num_elem);
                case 'reduce_cols_R_symmetric'
                    [out_min, out_r] = feval(kernel, out_min, out_r, M, b, nv, nv, num_elem);       
                otherwise
                    errror('Unkown kernel entry point')
            end
            out_r = out_r + 1;
            [min_val, j] = min(out_min);
            i = out_r(j);
        else
            error('?')
        end
            
        
        if abs(min_val(1)) < opt.Tol
            disp('Minimization complete because optimality tolerance was reached.')
            disp(['m : ', num2str(min_val)])
            break
        end
        
        x(i) = x(i) + 1;
        x(j) = x(j) - 1;
        b = b + 2*M(:,i) - 2*M(:,j);
    end
    
    if iter == opt.Iterations
        fprintf('Optimization terminated because the iterations budget %d was reached\n', opt.Iterations);
        fprintf('min_val = %g\n', min_val);
    end
  
end
    
function [x] = gridsearch2(M1, x, x0, M2, A2, y, y0, opt)
    % Minimize (A1*x - x0)^T M1 (A1*x - x0) + (A2*y - y0)^T M2 (A2*y - y0)

    y0_tilde = y + y0;
    if opt.UseGPU
        M1_tilde = gpuArray(single(M1));
        b1 = gpuArray(single(2*M1*(x - x0)));
        dm1 = gpuArray(single(diag(M1_tilde)));
        M2_tilde = gpuArray(single(A2'*M2*A2));
        b2 = gpuArray(single(2*A2'*M2'*(A2*x - y0_tilde)));
        dm2 = gpuArray(single(diag(M2_tilde)));
    else
        M1_tilde = M1;
        b1 = 2*M1*(x - x0);
        dm1 = diag(M1_tilde);
        M2_tilde = A2'*M2*A2;
        b2 = 2*A2'*M2'*(A2*x - y0_tilde);
        dm2 = diag(M2_tilde);
    end
    if opt.Kernel
        warning('gridsearch2 does not support kernel. Using bsx instead...');
    end
    for iter = 1:opt.Iterations    
        [min_val, i] = min(bsxfun(@plus, dm1+b1, (dm1-b1)') - 2*M1_tilde + ...
                           bsxfun(@plus, dm2+b2, (dm2-b2)') - 2*M2_tilde);
        [min_val, j] = min(min_val);
        i = i(j);
        
        if abs(min_val(1)) < opt.Tol
            disp('Minimization complete because optimality tolerance was reached.')
            disp(['m : ', num2str(min_val)])
            break
        end
        
        x(i) = x(i) + 1;
        x(j) = x(j) - 1;
        b1 = b1 + 2*M1_tilde(:,i) - 2*M1_tilde(:,j);
        b2 = b2 + 2*M2_tilde(:,i) - 2*M2_tilde(:,j);
        
        if opt.Histories
            %E_hist1(iter) = gather((A1*x-x0)'*M1*(A1*x-x0));
            %E_hist2(iter) = gather((A2*x-y0_tilde)'*M2*(A2*x-y0_tilde));
            %m_hist(iter) = gather(min_val);
            update_histories(iter, M1, x, x0, M2, A2, y, y0, min_val);
        end
    end
    
    if iter == opt.Iterations
        disp('Minimization complete because the max number of iterations was reached.')
        disp(['m : ', num2str(min_val)])
    end
    
    %if opt.Plot
    %    E_hist1 = E_hist1(1:iter-1);
    %    E_hist2 = E_hist2(1:iter-1);
    %    m_hist = m_hist(1:iter-1);
    %end
end
    
function [x] = gridsearch2_optimized(M1, x, x0, M2, A2, y, y0, opt)
    % Minimize (x-x0)^T M1 (x-x0) + (A2 x - y0_tilde)^T M2 (A2 x -
    % y0_tilde)
    
    y0_tilde = y + y0;
    if opt.UseGPU
        F = gpuArray(single(M1 + A2'*M2*A2));
        df = gpuArray(single(diag(F)));
        f = M1*x0 + A2'*M2*y0_tilde;
        g = gpuArray(single(2*F'*x - 2*f));
    else
        F = M1 + A2'*M2*A2;
        df = diag(F);
        f = M1*x0 + A2'*M2*y0_tilde;
        g = 2*F'*x - 2*f;
    end
    for iter = 1:opt.Iterations
        if (opt.UseGPU && opt.bsx) || ~opt.UseGPU
            [min_val, i] = min(bsxfun(@plus, df+g, (df-g)') - 2*F);
            [min_val, j] = min(min_val);
            i = i(j);
        elseif opt.Kernel
            [out_min, out_r] = feval(kernel, out_min, out_r, F, df, g, nv, nv, num_elem);
            out_r = out_r + 1;
            [min_val, j] = min(out_min);
            i = out_r(j);
        else
            error('?')
        end
        
        if abs(min_val(1)) < opt.Tol
            disp('Minimization complete because optimality tolerance was reached.')
            disp(['m : ', num2str(min_val)])
            break
        end
        
        x(i) = x(i) + 1;
        x(j) = x(j) - 1;
        g = g + 2*F(:, i) - 2*F(:, j);
        
        if opt.Histories
            %E_hist1(iter) = gather((A1*x-x0)'*M1*(A1*x-x0));
            %E_hist2(iter) = gather((A2*x-y0_tilde)'*M2*(A2*x-y0_tilde));
            %m_hist(iter) = gather(min_val);
            update_histories(iter, M1, x, x0, M2, A2, y, y0, min_val)
        end
    end
    
    if iter == opt.Iterations
        disp('Minimization complete because the max number of iterations was reached.')
        disp(['m : ', num2str(min_val)])
    end
end
    
function [x, y] = gridsearch3_optimized(M1, x, x0, M2, A2, y, y0, opt)
    % Minimize (x-x0)^T M1 (x-x0) + (A2 x - y0_tilde)^T M2 (A2 x -
    % y0_tilde)
    
    y0_tilde = y + y0;
    if opt.UseGPU
        F = gpuArray(single(M1 + A2'*M2*A2));
        df = gpuArray(single(diag(F)));
        f = M1*x0 + A2'*M2*y0_tilde;
        g = gpuArray(single(2*F'*x - 2*f));
    else
        F = M1 + A2'*M2*A2;
        df = diag(F);
        f = M1*x0 + A2'*M2*y0_tilde;
        g = 2*F'*x - 2*f;
    end
    for iter = 1:opt.Iterations
        if (opt.UseGPU && opt.bsx) || ~opt.UseGPU
            [min_val, i] = min(bsxfun(@plus, df+g, (df-g)') - 2*F);
            [min_val, j] = min(min_val);
            i = i(j);
        elseif opt.Kernel
            [out_min, out_r] = feval(kernel, out_min, out_r, F, df, g, nv, nv, num_elem);
            out_r = out_r + 1;
            [min_val, j] = min(out_min);
            i = out_r(j);
        else
            error('?')
        end
        
        if abs(min_val(1)) < opt.Tol
            disp('Minimization complete because optimality tolerance was reached.')
            disp(['m : ', num2str(min_val)])
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
            %E_hist1(iter) = gather((A1*x-x0)'*M1*(A1*x-x0));
            %E_hist2(iter) = gather((A2*x-y0_tilde)'*M2*(A2*x-y0_tilde));
            %k = [x; y];
            %res = TCODS(mesh, 'k', k, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', false);
            %E_outer_hist(iter) = res.miq_energy;
            %m_hist(iter) = gather(min_val);
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

function [Ztilde, Rtilde] = resistance_distance()
    eps = opt.JLEps;
    JLFac = opt.JLFac;
    k = round(JLFac * log(nv) / eps^2);
    Y = 2*(rand(k, ne) > 0.5) - 1;
    Y = (1/sqrt(k)) * Y;
    Y = Y * d0;
    Ztilde = (F \ Y')';
    if opt.UseGPU
        Ztilde = gpuArray(single(Ztilde));
        Rtilde = pdist(Ztilde', 'squaredeuclidean')'; % note that Rtilde is only the upper triangle part (1d vector)
        Ztilde = gather(Ztilde);
    else
        Rtilde = pdist(Ztilde', 'squaredeuclidean')';
    end
%     if opt.Colamd
%         Lchol = ichol(Lperm,struct('type','ict','droptol',1e-4,'michol','off'));
%     else
%         Lchol = ichol(L,struct('type','ict','droptol',1e-4,'michol','off'));
%     end
% 
%     eps = opt.JLEps;
%     JLFac = opt.JLFac;
%     k = round(JLFac * log(nv) / eps^2);
%     %k = round(JLFac * log(ne / eps^2));
%     % Approximate Z = Q d0 Lp
%     Q = 2*(rand(k, ne) > 0.5) - 1;
%     Q = (1/sqrt(k)) * Q;
%     Y = Q*d0;
%     if opt.UseGPU && opt.CholGPU
%         %Lchol = gpuArray(single(full(Lchol)));
%         %Ztilde = single(zeros(k, nv, 'gpuArray'));
%         Lchol = gpuArray(full(Lchol));
%         Ztilde = zeros(k, nv, 'gpuArray');
%         Y = gpuArray(Y);
%     else
%         Ztilde = zeros(k, nv);
%     end
%     if opt.Colamd
%         Y = Y(:, p);
%     end
%     if opt.ZtildeOneBS
%         Ztilde = (Lchol' \ (Lchol \ Y'))';
%     else
%         for ii = 1:k
%             yi = Y(ii, :);
%             [zi, flag] = pcg(Lperm, yi', 1e-6, 100, Lchol, Lchol');
%             %zi = Lchol' \ (Lchol \ yi');
%             %if opt.Colamd
%             %    zi = zi(pt);
%             %end
%             Ztilde(ii, :) = zi';
%         end
%     end
% 
%     if opt.Colamd
%         Ztilde = Ztilde(:, pt);
%     end
% 
%     if opt.UseGPU && opt.CholGPU
%         clear Y;
%     end
% 
%     if opt.SampleResistance
%         Rtilde = sparse(nv, nv);
%     else
%         if opt.UseGPU
%             Ztilde = gpuArray(single(Ztilde));
%             Rtilde = pdist(Ztilde', 'squaredeuclidean')'; % note that Rtilde is only the upper triangle part (1d vector)
%             Ztilde = gather(Ztilde);
%         else
%             Rtilde = pdist(Ztilde', 'squaredeuclidean')';
%         end
%     end
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