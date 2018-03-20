function [res] = ioq_wrapper(fp, varargin)
    global LOG;
    global OUT_FOLDER_LP;
    m = Mesh(fp);
    [~, meshname, ~] = fileparts(fp);
    nv = m.nV;
    ne = m.nE;
    log_and_print(LOG, 'fp : %s\r\n', fp);
    log_and_print(LOG, 'nv : %d\r\n', nv);
    GD = gpuDevice();
    
    p = inputParser;

	addOptional(p, 'Iterations', 1000);
    addOptional(p, 'NSingularities', []);
    addOptional(p, 'Laplacian', 'conn');
    addOptional(p, 'LaplacianPInv', []);
    addOptional(p, 'Plot', false);
    addOptional(p, 'Tol', 1e-6);
    addOptional(p, 'highg_method', []);
    addOptional(p, 'beta_P', []);
    addOptional(p, 'n_alternating', 1);
    addOptional(p, 'BlockSize', 20000);
    addOptional(p, 'UseGPU', true);
    addOptional(p, 'bsx', false);
    addOptional(p, 'kernel', true);
    addOptional(p, 'Verbose', true);
    addOptional(p, 'Histories', false);
    addOptional(p, 'Debug', false);
    addOptional(p, 'InvMethod', 'GPUInv');
    addOptional(p, 'JLEps', 'GPUInv');
    
    parse(p, varargin{:});
    opt = p.Results;
    JL_approx = strcmpi(opt.InvMethod, 'ApproxResistance');
    
    %if m.genus == 0 && ~strcmpi(opt.highg_method, 'genus0') && ~JL_approx
    %    error('skipping because mesh is genus 0 and method is for high genus only')
    %end
    
    switch opt.Laplacian
        case 'conn'
            Lp_path = fullfile(OUT_FOLDER_LP, [meshname, '_Lp_conn.mat']);
        case 'cot'
            Lp_path = fullfile(OUT_FOLDER_LP, [meshname, '_Lp_cot.mat']);
        otherwise
            error('Unkown laplacian type')
    end
    
    elapsed_inv = 0;
    % If Lp doesn't exist on disk, recompute it
    if ~exist(Lp_path, 'file') && ~JL_approx 
        fprintf('Creating L...\r\n');
        if strcmpi(opt.Laplacian, 'conn')
 			[d0, ~] = get_exterior_derivatives(m);
   			L = d0'*d0;
   		elseif strcmpi(opt.Laplacian, 'cot')
   			L = -cotmatrix(m.V, m.F);
        end
        
        try
            log_and_print(LOG, 'Trying gpu inv...\r\n');
            tic
            Lp = inv(gpuArray(single(L+1/nv))) - 1/nv;
            Lp = double(gather(Lp));
            if ~isempty(GD), wait(GD); end; elapsed_inv = toc
            alg_used = 'ioq_gpu_inv';
        catch ME
            fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
            log_and_print(LOG, 'GPUInv failed, trying block_inv_gpu...\r\n')
            try
                tic
                Lp = block_inv_gpu(full(L+1/nv), opt.BlockSize, true) - 1/nv;
                Lp = double(Lp);
                if ~isempty(GD), wait(GD); end; elapsed_inv = toc
                alg_used = 'ioq_block_inv';
            catch ME
                fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
                log_and_print(LOG, 'block_inv_gpu failed, trying block_inv_gpu without block mult\r\n');
                try
                    Lp = block_inv_gpu(full(L+1/nv), opt.BlockSize, false) - 1/nv;
                    Lp = double(Lp);
                    if ~isempty(GD), wait(GD); end; elapsed_inv = toc
                    alg_used = 'ioq_block_inv';
                catch ME
                    log_and_print(LOG, 'block_inv_gpu failed, trying invChol_mex...\r\n')
                    fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
                    tic
                    Lp = invChol_mex(full(L + 1/nv)) - 1/nv;
                    elapsed_inv = toc
                    alg_used = 'ioq_cpu_inv';
                end
            end
        end
        
        lap_type = opt.Laplacian;
        save(Lp_path, 'alg_used', 'elapsed_inv', 'lap_type', 'fp', 'Lp', '-v7.3')
    elseif ~JL_approx % load mesh from disk
        log_and_print(LOG, 'Loading Lp from %s...', Lp_path);
        tic
        load(Lp_path, 'alg_used', 'elapsed_inv', 'lap_type', 'fp', 'Lp');
        toc
        Lp = full(Lp);
    else
        alg_used = 'jl_ioq';
        Lp = [];
    end
    
    try
        log_and_print(LOG, 'Trying gpu IOQ...\r\n');
        [alpha_p, beta_p, elapsed_ioq] = IOQ_highgenus_gpu(m.V, m.F, ...
             'UseGPU', true, ...
             'LaplacianPInv', Lp, ...
             'Mesh', m, ...
             varargin{:});
         alg_used = [alg_used '_iter_gpu'];
    catch ME
        fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
        log_and_print(LOG, 'IOQ failed with GPU, trying without...\r\n');
        [alpha_p, beta_p, elapsed_ioq] = IOQ_highgenus_gpu(m.V, m.F, ...
             'UseGPU', false, ...
             'LaplacianPInv', Lp, ...
             'Mesh', m, ...
             varargin{:});
         alg_used = [alg_used '_iter_cpu'];
    end
    
    
%     try
%         log_and_print(LOG, 'Trying gpu inv...')
%         [alpha_p, beta_p, elapsed_total] = IOQ_highgenus_gpu(m.V, m.F, ...
%             'InvMethod', 'GPUInv', ...
%             varargin{:});
%         res.alg_used = 'ioq_gpu_inv_iter_gpu';
%     catch ME
%         log_and_print(LOG, '\tStack: %s\r\n', getReport(ME, 'extended'));
%         try 
%             log_and_print(LOG, 'GPUInv failed, trying block_inv_gpu...')
%             [alpha_p, beta_p, elapsed_total] = IOQ_highgenus_gpu(m.V, m.F, ...
%                 'InvMethod', 'GPUBlockInv', ...
%                 varargin{:});
%             res.alg_used = 'ioq_block_inv_iter_gpu';
%         catch
%             try
%                 log_and_print(LOG, 'block_inv_gpu failed, trying invChol_mex...')
%                 [alpha_p, beta_p, elapsed_total] = IOQ_highgenus_gpu(m.V, m.F, ...
%                     'InvMethod', 'CholMexInv', ...
%                     varargin{:});
%                 res.alg_used = 'ioq_cpu_inv_iter_gpu';
%             catch
%                 log_and_print(LOG, 'invChol_mex failed, trying without gpu...')
%                 [alpha_p, beta_p, elapsed_total] = IOQ_highgenus_gpu(m.V, m.F, ...
%                     'InvMethod', 'CholMexInv', ...
%                     'UseGPU', false, ...
%                     varargin{:});
%                 res.alg_used = 'ioq_cpu_inv_iter_cpu';
%             end
%         end            
%     end
    
    k = [alpha_p; beta_p];
    res.m = TCODS(m, 'k', k, 'f0', 1, 'theta0', 0, 'degree', 4, 'CreateFField', true, 'Duplicate', false, 'gConstraintVec', [1,0,0], 'Verbose', false);
    res.elapsed_total = elapsed_inv + elapsed_ioq;
    res.alg_used = alg_used;
    fprintf('Time : %g\n', elapsed_inv);
end

