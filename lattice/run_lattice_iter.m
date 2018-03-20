function res = run_lattice_iter(mesh, f0, theta0, degree, n_iter, opt)
    % Iterative {+1,-1} improvement heuristic.
    % See notes/2017-08-31 Iterative impr[3852].pdf
    %
    % Input:
    %   mesh - the mesh on which to calculate the direction field
    %   f0 - constrained faces
    %   theta0 - constraint angles
    %   degree - field degree
    %   n_iter - number of iterations
    %   opt - struct with the following boolean fields:
    %           graph - use graph Laplacian
    %           cot - use cot Laplacian
    %           gpu - use GPU
    %           chol - use cholesky decomposition
    %           eigs - use eigen decomposition
    %
    % Output:
    %   res - struct with the following fields:
    %           elapsed_iter - our heuristic elapsed time (in seconds)
    %           elapsed_TC - TC elapsed time
    %           elapsed - total elapsed time
    %           see TCODS for the rest.
    
    if ischar(mesh)
        fp = mesh;
        mesh = Mesh();
        mesh.loadTM(fp);
    end
    
    if ~isfield(opt, 'plot') , opt.plot = false         ; end
    if ~isfield(opt, 'inv_gpu')  , opt.inv_gpu = false  ; end
    if ~isfield(opt, 'iter_gpu') , opt.iter_gpu = false ; end
    if isfield(opt, 'graph')
        opt.cot = ~opt.graph;
    elseif isfield(opt, 'cot')
        opt.graph = ~opt.cot;
    else
        error('Lap type not was not given')
    end
    if ~isfield(opt, 'CUDAKernel') , opt.CUDAKernel = false ; end
    if ~isfield(opt, 'CUDAKernel_reduce_cols')
        opt.CUDAKernel_reduce_cols = false; 
    else
        opt.CUDAKernel = false;
    end
    if ~isfield(opt, 'TCODS'), opt.TCODS = false; end
    if ~isfield(opt, 'debug_CUDA'), opt.debug_CUDA = false; end
    if ~isfield(opt, 'tol'), opt.tol = 1e-5; end
    if ~isfield(opt, 'block_inv_gpu'), opt.block_inv_gpu = false; end
    if ~isfield(opt, 'lb'), opt.lb = 8000; end
    if ~isfield(opt, 'gConstraintVec'), opt.gConstraintVec = []; end

    nV = mesh.nV;

    if isfield(opt, 'alg_selector')
        opt.inv_gpu = opt.alg_selector.is_inv_gpu(nV);
        opt.block_inv_gpu = opt.alg_selector.is_inv_block(nV);
        opt.iter_gpu = opt.alg_selector.is_iter_gpu(nV);
        if ~opt.iter_gpu
            opt.CUDAKernel = false;
            opt.CUDAKernel_reduce_cols = false;
        end
    end

    global EPS;
    global LOG;
    if ~exist('LOG', 'var')
        LOG = -1;
    end
    if ~exist('EPS', 'var') || isempty(EPS)
        EPS = 1e-10;
    end
    
    using_gpu = opt.inv_gpu || opt.iter_gpu || opt.CUDAKernel || opt.CUDAKernel_reduce_cols || opt.block_inv_gpu;
    opt.iter_gpu = opt.iter_gpu || opt.CUDAKernel || opt.CUDAKernel_reduce_cols;
    
    if using_gpu
        try
            gd = gpuDevice();
        catch ME
            warning('gpuDevice() failed. Check that a GPU is available')
        end
    end
    
    
    
    log_and_print(LOG, '\t\tCreating L...\r\n');
    tic
    if opt.graph
        log_and_print(LOG, '\t\tUsing connectivity L...\r\n');
        [d0, ~] = get_exterior_derivatives(mesh);
        L = d0'*d0;
    elseif opt.cot
        log_and_print(LOG, '\t\tUsing cot L...\r\n');
        %L = lap(mesh);
        %L = -diag(1 ./ mesh.AV)*cotmatrix(mesh.V, mesh.F);
        %L = diag(1 ./ mesh.AV)*lap_cot(mesh);
        %options.symmetrize = 1;
        %options.normalize = 1;
        %L = diag(1 ./ (mesh.AV ./ norm(mesh.AV))) * compute_mesh_laplacian(mesh.V', mesh.F', 'conformal', options);
        %L = compute_mesh_laplacian(mesh.V', mesh.F', 'conformal', options);
        
        %gradop = grad(mesh);
        %divop = div(mesh, gradop);
        %L = -divop * gradop;
        
        L = -cotmatrix(mesh.V, mesh.F);
        %[~, L] = grad(mesh);
    end
    L = single(full(L));
    
    if exist('gd', 'var')
        wait(gd);
    end
    elapsed_lap = toc();
    
    %xi = 2 - 2*mesh.genus;
     
    log_and_print(LOG, '\t\tInverting L...\r\n');
    tic
    if opt.inv_gpu
        if opt.cot
            Lp = inv(gpuArray(L));
            %Lp = 0.25 * Lp * mesh.Gf * Lp' * mesh.Gv;
            %Lp = 0.25*mesh.Gf * (gpuArray(full(L')) \ full(mesh.Gv));
        else
            Lp = inv(gpuArray(L + 1/nV)) - 1/nV;
        end
    elseif opt.block_inv_gpu
        if opt.cot
            Lp = 0.25*mesh.Gv * block_inv_gpu(L)' * mesh.Gv;
        else
            Lp = block_inv_gpu(L + 1/nV, opt.lb);
        end
        try
            Lp = gpuArray(Lp);
        catch ME
            log_and_print(LOG, '\t\tLp of size %d does not fit on GPU. Using cpu loop instead.\r\n', size(Lp,1));
            opt.iter_gpu = false;
        end
    else
        Lp = invChol_mex(L + 1/nV) - 1/nV;
        if opt.iter_gpu
            Lp = gpuArray(Lp);
        end
    end
    
    if exist('gd', 'var')
        wait(gd); 
    end
    elapsed_inv = toc();
    
    % Rather disappointing
    %Lp = qrginv(L);

    log_and_print(LOG, '\t\tSetup...\r\n');
    
    tic
    Ad = get_gaussian_curvature(mesh);
    k0 = (2/pi)*Ad;
    if ~isfield(opt, 'n_singularities')
        c = round(abs(sum(k0))); % = 4*xi
    else
        c = opt.n_singularities;
    end
    
    k = zeros(nV, 1);
    n_pos_sing = (c + round(sum(k0))) / 2;
    n_neg_sing = (c - round(sum(k0))) / 2;
    inds_pos = randperm(nV, n_pos_sing);
    inds_neg = randperm(nV, n_neg_sing);
    k(inds_pos) = 1;
    k(inds_neg) = -1;
    
    dlp = diag(Lp);
    b = 2*(Lp*(k-k0));
    
    if opt.iter_gpu && opt.CUDAKernel
        % Setup the kernel
        kernel = parallel.gpu.CUDAKernel('lattice_inner_loop.ptx', 'lattice_inner_loop.cu', 'lattice_inner_loop');
        num_elem = nV^2;
        kernel.ThreadBlockSize = [kernel.MaxThreadsPerBlock, 1, 1];
        kernel.GridSize = [ceil(nV / kernel.MaxThreadsPerBlock), nV];       
        out = zeros(size(Lp), 'like', Lp);
        if ~opt.inv_gpu
            dlp = gpuArray(dlp);
            b = gpuArray(b);
        end
    elseif opt.iter_gpu && opt.CUDAKernel_reduce_cols
        % Setup the kernel
        kernel = parallel.gpu.CUDAKernel('lattice_inner_loop.ptx', 'lattice_inner_loop.cu', 'reduce_cols');
        num_elem = nV^2;
        kernel.ThreadBlockSize = [1, kernel.MaxThreadsPerBlock, 1];
        kernel.GridSize = [1, ceil(nV / kernel.MaxThreadsPerBlock)];
        out_min = single(zeros(1, nV, 'gpuArray'));
        out_r = uint32(zeros(1, nV, 'gpuArray'));
        if ~opt.inv_gpu
            dlp = gpuArray(dlp);
            b = gpuArray(b);
        end
    elseif opt.inv_gpu && ~opt.iter_gpu
        Lp = gather(Lp);
        dlp = gather(dlp);
        b = gather(b);
    end

    if exist('gd', 'var')
        wait(gd); 
    end
    elapsed_setup = toc();
    
    assert(n_pos_sing + n_neg_sing == c)
    assert(abs(sum(k)-sum(k0)) < EPS)
    
    log_and_print(LOG, '\t\tLoop...\r\n');
    tic
    for iter = 1:n_iter
        if opt.iter_gpu && opt.CUDAKernel
            out = feval(kernel, out, Lp, dlp, b, nV, nV, num_elem);
            [m, i] = min(out(:));
            [i, j] = ind2sub(size(out), i);
        elseif opt.iter_gpu && opt.CUDAKernel_reduce_cols
            [out_min, out_r] = feval(kernel, out_min, out_r, Lp, dlp, b, nV, nV, num_elem);
            out_r = out_r + 1;
            [m, j] = min(out_min);
            i = out_r(j);
        else
            [m, i] = min(bsxfun(@plus, dlp+b, (dlp-b)') - 2*Lp);
            [m, j] = min(m);
            i = i(j);
        end
        
        if opt.inv_gpu && opt.debug_CUDA
            %out2 = bsxfun(@plus, dlp+b, (dlp-b)') - 2*Lp;
            %[m2, i2] = min(out2(:));
            %[i2, j2] = ind2sub(size(out2), i2);
            %check_norm('out', 'out2');
            %check_norm('m', 'm2');
            
            if opt.CUDAKernel_reduce_cols
                out2 = bsxfun(@plus, dlp+b, (dlp-b)') - 2*Lp;
                [out_min2, out_r2] = min(out2);
                [m2, j2] = min(out_min2);
                i2 = out_r2(j2);
                check_norm('out_min', 'out_min2');
                check_norm('double(out_r)', 'double(out_r2)');
                check_norm('m', 'm2');
                check_norm('double(i)', 'i2');
                check_norm('double(j)', 'j2');
            end
        end
        
        if abs(m(1)) < opt.tol
            break;
        end

        k(i) = k(i) + 1;
        k(j) = k(j) - 1;
        
        b = b + 2*Lp(:,i) - 2*Lp(:,j);
        
        if opt.plot
            E_lat_hist(iter) = gather((k-k0)'*Lp*(k-k0));
            m_hist(iter) = gather(m);
            %min_D_hist(iter) = gather(min(D(:)));
            %min_R_hist(iter) = gather(min(R(:)));
            if opt.TCODS
                inds = find(abs(k) > EPS);
                S = [inds, k(inds) / degree];
                tmp_res = TCODS(mesh, S, f0, theta0, degree, false);
                E_MIQ_hist(iter) = tmp_res.E;
            end
        end
    end
    
    if exist('gd', 'var')
        wait(gd); 
    end
    elapsed_iter = toc;
    
    if opt.plot
        figure
        xx = 1:length(E_lat_hist);
        plot(xx, E_lat_hist, '-r', ...
            xx, m_hist, '-b', ...
            xx, E_MIQ_hist, 'm')
            %xx, min_D_hist, '-m', ...
            %xx, min_R_hist, '-g')
        legend('E lat', ...
            'm', ...
            'E MIQ')
            %'min(D)', ...
            %'min(R)')
    end
    
%     tic
    if opt.TCODS
        log_and_print(LOG, '\t\tTCODS...\r\n');
        inds = find(abs(k) > EPS);
        S = [inds, k(inds) / degree];
        tic
        % TODO mex crashes matlab for some reason
        res = TCODS(mesh, S, f0, theta0, degree, ...
            'Mex', false, ...
            'gConstraintVec', opt.gConstraintVec);
        res.elapsed_TCODS = toc();
        res.S = S;
        log_and_print(LOG, '\t\tE: %g\r\n', res.E);
    end
     
    res.k = k;
    res.elapsed_setup_lap = elapsed_lap;
    res.elapsed_setup = elapsed_setup;
    res.elapsed_inv = elapsed_inv;
    res.elapsed_iter = elapsed_iter;
    res.elapsed_total = elapsed_lap + elapsed_setup + elapsed_inv + elapsed_iter;
    res.nV = mesh.nV;
    res.nF = mesh.nF;

    % Algorithm used
    if opt.iter_gpu
        if opt.inv_gpu
            res.alg_used = 'ioq_gpu_inv_iter_gpu';
        elseif opt.block_inv_gpu
            res.alg_used = 'ioq_block_inv_iter_gpu';
        else
            res.alg_used = 'ioq_cpu_inv_iter_gpu';
        end
    else
        if opt.inv_gpu
            error('Should not reach this')
        elseif opt.block_inv_gpu
            res.alg_used = 'ioq_block_inv_iter_cpu';
        else
            res.alg_used = 'ioq_cpu_inv_iter_cpu';
        end
    end
     
%     res.elapsed_TC = toc;
%     res.elapsed_iter = elapsed_iter;
%     
%     log_and_print(LOG, 'E = %g\r\n', res.E);
%     log_and_print(LOG, 'Done.\r\n');
    
%     res.elapsed = res.elapsed_iter + res.elapsed_TC;
    res.title = {'Iterative', ...
        sprintf('Total time: %g', res.elapsed_total)};
    if opt.TCODS
        res.title = {res.title{:}, sprintf('$E_{MIQ} = %g$', res.E)};
    end