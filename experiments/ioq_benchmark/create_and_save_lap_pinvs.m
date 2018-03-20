% Invert 

ioq_benchmark_setup;
filepaths = get_filepaths(data_folder, data_ext);
n_files = numel(filepaths);
GD = gpuDevice();
BLOCK_SIZE = 20000;

print_header(sprintf('base folder: %s', base_folder));

for i = 1:n_files
    fp = filepaths{i};
    print_header(fp);
    m = Mesh(fp); nv = m.nV;
    [~, meshname, ~] = fileparts(fp);
    fprintf('fp : %s\n', fp);
    fprintf('nv : %d\n', nv);
    
    Lp_path_conn = fullfile(OUT_FOLDER_LP, [meshname, '_Lp_conn.mat']);
    Lp_path_cot = fullfile(OUT_FOLDER_LP, [meshname, '_Lp_cot.mat']);
    if exist(Lp_path_conn, 'file') && exist(Lp_path_cot, 'file')
        disp('Skipping...')
        continue
    end
    
    % Calc and save conn pinv
    fprintf('Inverting conn lap...\n');
    d0 = get_exterior_derivatives(m);
    L = d0' * d0;
    
    try
        disp('Trying gpu inv...')
        tic
        Lp = inv(gpuArray(single(L+1/nv))) - 1/nv;
        Lp = double(gather(Lp));
        if ~isempty(GD), wait(GD); end; elapsed_inv = toc
        alg_used = 'ioq_gpu_inv';
    catch ME
        disp('GPUInv failed, trying block_inv_gpu...')
        fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
        try
            tic
            Lp = block_inv_gpu(full(L+1/nv), BLOCK_SIZE, true) - 1/nv;
            Lp = double(Lp);
            if ~isempty(GD), wait(GD); end; elapsed_inv = toc
            alg_used = 'ioq_block_inv';
        catch ME
            disp('block_inv_gpu failed, trying block_inv_gpu without block mult');
            fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
            try
                Lp = block_inv_gpu(full(L+1/nv), BLOCK_SIZE, false) - 1/nv;
                Lp = double(Lp);
                if ~isempty(GD), wait(GD); end; elapsed_inv = toc
                alg_used = 'ioq_block_inv';
            catch ME
                disp('block_inv_gpu failed, trying invChol_mex...')
                fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
                try
                    tic
                    Lp = invChol_mex(full(L + 1/nv)) - 1/nv;
                    elapsed_inv = toc
                    alg_used = 'ioq_cpu_inv';
                catch
                    disp(['Could not calculate Lp for ', fp])
                    fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
                    continue
                end
            end
        end
    end
    
    inds = find(Lp<1e-10);
    fprintf('Sparsifying Lp by %g%%\n', 100*(1 - length(inds)/nv^2));
    Lp(inds) = 0;
    Lp = sparse(Lp);
    lap_type = 'conn';
    Lp_path = fullfile(OUT_FOLDER_LP, [meshname, '_Lp_conn.mat']);
    save(Lp_path, 'alg_used', 'elapsed_inv', 'lap_type', 'fp', 'Lp', '-v7.3')
    
    % Calc and save cot pinv
    fprintf('Inverting cot lap...\n');
    L = -cotmatrix(m.V, m.F);
    
    try
        disp('Trying gpu inv...')
        tic
        Lp = inv(gpuArray(single(L+1/nv))) - 1/nv;
        Lp = double(gather(Lp));
        if ~isempty(GD), wait(GD); end; elapsed_inv = toc
        alg_used = 'ioq_gpu_inv';
    catch ME
        disp('GPUInv failed, trying block_inv_gpu...')
        fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
        try
   
            tic
            Lp = block_inv_gpu(full(L+1/nv), BLOCK_SIZE, true) - 1/nv;
            Lp = double(Lp);
            if ~isempty(GD), wait(GD); end; elapsed_inv = toc
            alg_used = 'ioq_block_inv';
        catch ME
            disp('block_inv_gpu failed, trying block_inv_gpu without block mult...')
            fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
            try   
                Lp = block_inv_gpu(full(L+1/nv), BLOCK_SIZE, false) - 1/nv;
                Lp = double(Lp);
                if ~isempty(GD), wait(GD); end; elapsed_inv = toc
                alg_used = 'ioq_block_inv';
            catch ME
                disp('block_inv_gpu failed, trying invChol_mex...')
                fprintf('\tStack: %s\r\n', getReport(ME, 'extended'));
                try
                    tic
                    Lp = invChol_mex(full(L + 1/nv)) - 1/nv;
                    elapsed_inv = toc
                    alg_used = 'ioq_cpu_inv';
                catch
                    disp(['Could not calculate Lp for ', fp]);
                    continue
                end
            end
        end
    end
    
    inds = find(Lp<1e-10);
    fprintf('Sparsifying Lp by %g%%\n', 100*(1 - length(inds)/nv^2));
    Lp(inds) = 0;
    Lp = sparse(Lp);
    lap_type = 'cot';
    Lp_path = fullfile(OUT_FOLDER_LP, [meshname, '_Lp_cot.mat']);
    save(Lp_path, 'alg_used', 'elapsed_inv', 'lap_type', 'fp', 'Lp', '-v7.3')
end