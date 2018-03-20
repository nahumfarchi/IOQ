function [Rtilde, Ztilde] = resistance_distance(m, eps, JLFac, useGPU, k)
%[Ztilde, Rtilde] = resistance_distance(m, eps, JLFac, useGPU)
    if nargin < 2
        eps = 0.5;
    end
    if nargin < 3
        JLFac = 24;
    end
    if nargin < 4
        useGPU = false;
    end
    
    nv = m.nV; ne = m.nE;
    d0 = get_exterior_derivatives(m);
    L = d0' * d0;
    F = factorize(L);
    
    if eps == 0       
        %Lp = F \ eye(nv);

        %dl = diag(Lp);
        if useGPU
            %Lp = inv(gpuArray(single(full(L))+1/nv)) - 1/nv;
            Lp = inv(single(gpuArray(full(L)+1/nv))) - 1/nv;
            %Lp = gpuArray(single(full(Lp)));
            dl = gpuArray(single(full(diag(Lp))));
        else
            Lp = invChol_mex(single(full(L))+1/nv) - 1/nv;
            dl = diag(Lp);
        end
        Rtilde = bsxfun(@plus, dl, dl') - 2*Lp;
        Ztilde = [];
    else
        if ~exist('k', 'var')
            k = round(JLFac * log(nv) / eps^2);
        end
        if k > nv
            k = nv;
        end
        Y = 2*(rand(k, ne) > 0.5) - 1;
        Y = (1/sqrt(k)) * Y;
        %Y = randi([-1, 1], k, ne) / sqrt(k);
        Y = Y * d0;
        Ztilde = (F \ Y');
        % note that Rtilde is only the upper triangle part (1d vector)
        if useGPU
            tic
            Ztilde = single(gpuArray(Ztilde));
            %Ztilde = gpuArray(Ztilde);
            if nv < 60000
                Rtilde = pdist(Ztilde, 'squaredeuclidean')'; 
            else
                Rtilde = pdist_block();
            end
            %Ztilde = gather(Ztilde);
        else
            Rtilde = pdist(Ztilde, 'squaredeuclidean')';
        end
        
        %Rtilde2 = pdist_block();
        %assert(norm(squareform(Rtilde) - Rtilde2, 'fro') < 1e-4)
        %assert(norm(Rtilde - Rtilde2, 'fro') < 1e-10)
    end
    
function res = pdist_block()
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
    
end