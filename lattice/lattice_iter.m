function res = lattice_iter(mesh, f0, theta0, degree, n_iter, opt)
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
    
    if ~isfield(opt, 'plot') , opt.plot = false  ; end
    if ~isfield(opt, 'gpu')  , opt.gpu = false; opt.only_Lp_on_gpu=false ; end
    if ~isfield(opt, 'graph'), opt.graph = false ; opt.cot = true; end
    if ~isfield(opt, 'cot')  , opt.graph = true  ; opt.cot = false; end
    if ~isfield(opt, 'chol') , opt.chol = false  ; end
    if ~isfield(opt, 'eigs') , opt.eigs = false  ; end
       
    if opt.gpu || opt.only_Lp_on_gpu
        gpuDevice(1);
    end
    
    tic
    LOG = -1;
    EPS = 1e-10;
    
    nV = mesh.nV;
    
    log_and_print(LOG, 'Creating L...\r\n');
    if opt.graph
        [d0, ~] = get_exterior_derivatives(mesh);
        L = d0'*d0;
    elseif opt.cot
        %L = lap(mesh);
        L = lap_cot(mesh);
    end
    Ad = get_gaussian_curvature(mesh);
    %xi = 2 - 2*mesh.genus;
     
    % faster way to get the pseudo inverse   
    log_and_print(LOG, 'Inverting L...\r\n');
    if opt.chol
        Lp = invChol_mex(L + 1/nV) - 1/nV;
        if opt.gpu
            Lp = gpuArray(single(Lp));
        end
    elseif opt.eigs
        [V, Dp] = eigs(L, opt.n_eig, 'sm');
        Dp = diag(Dp);
        Dp(Dp < EPS) = inf;
        Dp = diag(1 ./ Dp);
        Lp = V*Dp*V';
        Lp = gpuArray(full(single(Lp)));
    elseif opt.gpu
        Lp = inv(gpuArray(single(L + 1/nV))) - 1/nV;
    elseif opt.only_Lp_on_gpu
        Lp = gpuArray(invChol_mex(single(L+1/nV)) - 1/nV);
    else
        Lp = inv(L + 1/nV) - 1/nV;
    end
    
    % Rather disappointing
    %Lp = qrginv(L);

    k0 = (2/pi)*Ad;
    %c = 4*xi;
    c = round(sum(k0));
    
    k = zeros(nV, 1);
    inds = randperm(nV, c);
    k(inds) = 1;
    
    % Resistance distance matrix
    log_and_print(LOG, 'Calculating resistance matrix...\r\n');
    
    if opt.gpu
        % Not clear which one is better
        %Ldii = repmat(diag(Lp), 1, nV);
        %R = gpuArray(Ldii + Ldii' - 2*Lp);
        %tic
        dlp = diag(Lp); 
        R = bsxfun(@plus,dlp,dlp') - 2*Lp; 
        %toc
        b = gpuArray(2*Lp*(k-k0)); % runs out of memory with backslash (on large meshes)
    elseif opt.only_Lp_on_gpu
        dlp = diag(Lp); 
        %R = gather(bsxfun(@plus,dlp,dlp') - 2*Lp); 
        b = 2*(Lp*(k-k0)); % runs out of memory with backslash (on large meshes)
    else
        %Ldii = repmat(diag(Lp), 1, nV);
        %R = gather(Ldii + Ldii' - 2*Lp);
        %b = 2*Lp*(k-k0);
        dlp = diag(Lp); 
        R = bsxfun(@plus,dlp,dlp') - 2*Lp; 
        b = 2*(Lp \ (k-k0));
    end
    
    log_and_print(LOG, 'Iterative improvement...\r\n');
    for iter = 1:n_iter
        if opt.only_Lp_on_gpu
            [m, i] = min(bsxfun(@plus, dlp, dlp') - 2*Lp + bsxfun(@minus, b, b'));
            [m, j] = min(m);
            i = i(j);
        else
            D = gather(bsxfun(@minus, b, b'));
    
            Z = D + R;
            [m, i] = min(Z(:));
        end
        
        if abs(m) < EPS
            break;
        end
        
        [i, j] = ind2sub(size(Z), i);
        k(i) = k(i) + 1;
        k(j) = k(j) - 1;
        
        b = b + 2*Lp(:,i) - 2*Lp(:,j);
        
        if opt.plot
            E_lat_hist(iter) = gather((k-k0)'*Lp*(k-k0));
            m_hist(iter) = gather(m);
            min_D_hist(iter) = gather(min(D(:)));
            min_R_hist(iter) = gather(min(R(:)));
        end
        
        if ~opt.only_Lp_on_gpu
            clear Z;
        end
    end
    wait(gpuDevice); elapsed_iter = toc;
    
    
    
    if opt.plot
        figure
        xx = 1:length(E_lat_hist);
        plot(xx, E_lat_hist, '-r', ...
            xx, m_hist, '-b', ...
            xx, min_D_hist, '-m', ...
            xx, min_R_hist, '-g')
        legend('E lat', ...
            'm', ...
            'min(D)', ...
            'min(R)')
    end
    
%     tic
     inds = find(abs(k) > EPS);
     S = [inds, k(inds) / degree];
     res = TCODS(mesh, S, f0, theta0, degree, false);
     
     res.S = S;
     res.k = k;
     res.elapsed_iter = elapsed_iter;
%     res.elapsed_TC = toc;
%     res.elapsed_iter = elapsed_iter;
%     
%     log_and_print(LOG, 'E = %g\r\n', res.E);
%     log_and_print(LOG, 'Done.\r\n');
    
%     res.elapsed = res.elapsed_iter + res.elapsed_TC;
    res.title = {'Iterative', ...
        sprintf('$E_{MIQ} = %g$', res.E), ...
        sprintf('Elapsed time: %g', res.elapsed_iter)};
    
    function chunk_inner_loop()
        gd = gpuDevice();
        n = size(Lp, 1);
        n_half = floor(n / 2);
        mem_required = n_half^2 * sizeof(single);
        mem_avail = gd.AvailableMemory;
        if mem_required < mem_avail
            X11 = bsxfun(dlp(1:n_half) + b(1:n_half), (dlp(1:n_half) - b(1:n_half))');
            X12 = bsxfun(dlp(1:n_half) + b(1:n_half), (dlp(n_half+1:end) - b(n_half+1:end))');
            X21 = bsxfun(dlp(n_half+1:end) + b(n_half+1:end), (dlp(1:n_half) - b(1:n_half))');
            X22 = bsxfun(dlp(n_half+1:end) + b(n_half+1:end), (dlp(n_half+1:end) - b(n_half+1:end))');
            
            [m11, i11] = min(X11(:));
            [i11, j11] = ind2sub(size(X11), i11);
            
            [m12, i12] = min(X12(:));
            [i12, j12] = ind2sub(size(X12), i12);
            j12 = n_half + j12;
            
            [m21, i21] = min(X21(:));
            [i21, j21] = ind2sub(size(X21), i21);
            i21 = n_half + i21;
            
            [m22, i22] = min(X22(:));
            [i22, j22] = ind2sub(size(X22), i22);
            i22 = n_half + i22;
            j22 = n_half + j22;
        else % Split again
            X11 = 
        end
    end
end