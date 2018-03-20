% Find the largest mesh that can be inverted on the gpu before running out
% of memory.

filepaths = get_filepaths('../data/bunnies', '.off');
gd = gpuArray();

best_nF = -inf;
for fp = filepaths
    disp(fp)
    m = Mesh();
    m.loadTM(fp{:})
    nV = m.nV;
    
    Lsp = -cotmatrix(m.V, m.F); % sparse
    Lfl = full(Lsp); % full
    Lsp_gpu = gpuArray(Lsp);
    Lfl_gpu = gpuArray(Lfl);

    try
        tic; inv(Lsp_gpu); wait(gd); toc
    catch ME
        fprintf('Largest: %d\n', best_nF);
        break;
    end
    best_nF = m.nF;
end

