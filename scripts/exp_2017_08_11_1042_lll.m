function [] = exp_2017_08_01_1234_globally_optimal()
ME = [];
try
%MESHES = {'sphere_s0.off'};
ABOUT = 'LLL Lattice reduction\r\n\r\n';
OUT_FOLDER_NAME = 'LLL';

%MESHES = {'sphere_s0.off', ...
%    'bumpy.off', ...
%    'round_cuber.off', ...
%    'rounded_cube_keenan.off', ...
%    'cow.off', ...
%    'bunny.off', ...
%    'bunny2.off', ...
%    'phands.off', ...
%    'torus_fat_r2.off'};

%MESHES = {'bunny2.off', 'bunny.off', 'rounded_cube_keenan.off', 'round_cuber.off'};
%MESHES = {'rounded_cube_keenan.off'};
%MESHES = {'round_cuber.off'};
%MESHES = {'sphere_s0.off'};
MESHES = {'sphere_s0.off'};

%VIEW_ANGLE = [];

VERBOSE = true;
EPS = 1e-9;

theta0 = 0;

RUN_LIBIGL_MIQ = true;
PLOT = false;
SAVE = false;

if PLOT && SAVE
    OUT_FOLDER = create_time_stamped_folder(fullfile('..', 'results'), ...
        OUT_FOLDER_NAME, ...
        true, ...
        false);
    LOG = fopen(fullfile(OUT_FOLDER, 'log.txt'), 'w');
else
    LOG = -1;
end

log_and_print(LOG, '%s\r\n', mfilename('fullpath'));
log_and_print(LOG, ABOUT);
log_and_print(LOG, 'logid : %d\r\n', LOG);

k = 1;
for fname = MESHES
    log_and_print(LOG, 'Loading %s...\r\n', fname{:});
    p = find_data_folder();
    fp = fullfile(p, fname{:});

    mesh = Mesh();
    mesh.loadTM(fp);
    
    [d0, d1] = get_exterior_derivatives(mesh);
    Ad = get_gaussian_curvature(mesh);
    degree = 4;
    f0 = [1]; % Starting face
    v0 = mesh.V(mesh.F(f0, 2), :) - mesh.V(mesh.F(f0, 1), :);
    xi = 2 - 2*mesh.genus;

    log_and_print(LOG, 'degree : %d\r\n', degree);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % libigl MIQ (greedy rounding) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if RUN_LIBIGL_MIQ
        log_and_print(LOG, '\r\nSolving with libigl MIQ...\r\n');
        MIQ = nrosy_wrapper(fp, f0, v0, degree);
        log_and_print(LOG, '%s\r\n', MIQ.result);
        log_and_print(LOG, 'Status: %d\r\n', MIQ.status);
        log_and_print(LOG, 'E: %g\r\n', MIQ.E);
    end

    % |Mk-t|^2
    N_EIGENVALUES = mesh.nV;
    W = d0'*d0;
    [V,D] = eig(full(W));
    [D, sort_eigen] = sort(diag(D));
    V = V(:, sort_eigen);
    D = diag(D);

    D = D(1:N_EIGENVALUES, 1:N_EIGENVALUES);
    V = V(:, 1:N_EIGENVALUES);

    Dinv = diag(D);
    Dinv(Dinv < EPS) = inf;
    Dinv = 1 ./ Dinv;
    Dinv = diag(Dinv);
    Winv = V*Dinv*V';

    % |M*x-t|
    M = V*sqrt(Dinv)*V';
    X = V*Dinv*V';
    b = (2/pi)*Ad;
    t = (2/pi)*M*Ad;
    %M = V*D*V';
    %t = V*D*sqrt(D)*t;
    m = 4*xi;
    
    %delta_ = 0.75;
    %tic
    %Ar = LLL_reduction(M(:,1:end-1), delta);
    %toc
    %
    %tic
    %k = babai(triu(Ar),t);
    %toc
    
    %k = LLL([M', ones(size(M',1),1)], [t', 4*xi]);
    %k = lattice_solve(M', t');
    %k = k(:);
    
%     % LLL_FP treats *rows* as basis vectors, while we treat *columns*.
%     M = LLL_FP(M')';
     %M = LLL_reduction(M);
     M_reduced = CLLL(M);
     M_B = M_reduced(:, 1);
     M_I = M_reduced(:, 2:end);
     c = 4*xi;
     A = M_I - M_B*ones(1, numel(M_B)-1);
     y = M_B*c - (M_reduced \ t);
     [Q, R] = qr(A);
     [Q2, R2] = qr(Q'*A); 
     yy = R2\y;
     k = babai(R2, yy);
%     k2 = babai2(A, y);
%     %[ksorted, inds] = sort(k);
%     %kmax = k(inds(end:-1:end-4));
%     %kmin = k(inds(1:4));
    
%     M_B = M(:, 1);
%     M_I = M(:, 2:end);
%     c = 4*xi;
%     A = M_I - M_B*ones(1, numel(M_B)-1);
%     y = M_B*c - t;
%     k = lattice_solve(A', y')';
%     k = [k; c-ones(1, numel(k))*k];
    
    %inds = find(abs(k) > EPS);
    %inds = inds(end:-1:end-7);
    %inds = inds(1:8);
    inds = find(abs(k) > EPS);
    S = [inds, k(inds)];
    LATTICE = TCODS(mesh, S, f0, theta0, degree, VERBOSE);
    fprintf('E = %g\n', LATTICE.E);
    
    LATTICE.title = {'Lattice lll', sprintf('$E_{MIQ} = %g$', LATTICE.E)};
    MIQ.title = {'MIQ ', sprintf('$E_{MIQ} = %g$', MIQ.E)};
    
    figure(1); clf(1)
    subplot(221)
    MeshVis.plot(mesh, ...
        'nrosy', LATTICE, ...
        'Title', LATTICE.title);
    subplot(222)
    MeshVis.plot(mesh, ...
        'nrosy', MIQ, ...
        'Title', MIQ.title);
end

% A hack since matlab does not have a finally clause
catch ME
end
% Close open resources
if LOG > 0
    fclose(LOG);
end
if ~isempty(ME)
    rethrow(ME);
end

    function A = LLL_reduction(A,delta)
        % LLL lattice reduction
        % Matlab implementation by K. Shum
        % The input A is a square matrix. Columns of input A are the basis vectors 
        % The output A is a square matrix, whose columns are LLL-redcued basis
        % vectors

        if nargin == 1
            delta = .75; % the default value of the parameter delta
        end

        m_ = length(A); % the dimension of the vector space
        B = zeros(m_,m_); % Columns of B are the vectors after the Gram-Schmidt process
        mu = zeros(m_,m_); % The matrix mu stores the GS coefficients
        M_ = zeros(1,m_); % M(i) is the norm squared of the i-th column of B

        % Gram-Schmidt orthogonalization
        B(:,1) = A(:,1); % Set the first column of B as the same
            % as the first column of A
        M_(1) = dot(B(:,1), B(:,1));
        for i = 2:m_
            mu(i,1:(i-1)) = (A(:,i)'* B(:,1:(i-1))) ./ M_(1:(i-1));
            B(:,i) = A(:,i) - B(:,1:(i-1))*mu(i,1:(i-1))';
            M_(i) = dot(B(:,i), B(:,i));
        end

        k_ = 2;
        while k_ <= m_

           % Size reduction
           for i = (k_-1):-1:1
               q = round(mu(k_,i));
               if q ~=0
                   A(:,k_) = A(:,k_) - q*A(:,i);  % size-reduce the k-th basis vector
                   mu(k_,1:i) =  mu(k_,1:i) - q*[mu(i,1:(i-1)) 1]; % update the GS coefficients
               end
           end

           % Check the Lovasz condition
           if dot(B(:,k_),B(:,k_)) >= (delta - abs(mu(k_,k_-1))^2)* dot(B(:,(k_-1)),B(:,(k_-1)))
               k_ = k_+1; % increment k_ if the Lovasz condition holds
           else
               % If the Lovasz condition fails,
               % swap the k_-th and (k_-1)-st basis vector
               v = A(:,k_); A(:,k_) = A(:,k_-1); A(:,k_-1) = v;

               % update the Gram-Schmidt coefficients
               for s = (k_-1):k_
                   mu(s,1:(s-1)) = (A(:,s)'* B(:,1:(s-1))) ./ M_(1:(s-1));
                   B(:,s) = A(:,s) - B(:,1:(s-1))*mu(s,1:(s-1))';
                   M_(s) = dot(B(:,s), B(:,s));
               end
               mu((k_+1):m_,(k_-1):k_) = (A(:,(k_+1):m_)'* B(:,(k_-1):k_)) / diag(M_((k_-1):k_));

               if k_ > 2
                  k_ = k_-1;
               end
            end
        end
    end

    function z_hat_ = babai(R_,y_)
    %%
    %   compute the Babai estimation
    %   find a sub-optimal solution for min_z ||R*z-y||_2
    %   R - an upper triangular real matrix of n-by-n
    %   y - a real vector of n-by-1
    %   z_hat - resulting integer vector
    %
        %n_=length(y_);
        n_ = size(R_, 2);
        z_hat_ = zeros(n_,1);
        z_hat_(n_) = round(y_(n_)./R_(n_,n_));
        for k_=n_-1:-1:1
            par_=R_(k_,k_+1:n_)*z_hat_(k_+1:n_);
            ck_=(y_(k_)-par_)./R_(k_,k_);
            z_hat_(k_)=round(ck_);
        end
    end

    function z_hat_ = babai2(B_, t_)
        n_ = size(B_, 2);
        R_ = orth(B_);
        b_ = t_ / norm(t_);
        for j_ = n_:-1:1
            bj_ = R_(:, j_);
            cj_ = round(dot(b_, bj_) / dot(bj_, bj_));
            b_ = b_ - cj_ * bj_;
        end
        z_hat_ = t_ - b_;
    end

    function B_reduced = CLLL(B_)

        % LLL algorithm of lattice reduction for complex-valued lattices

        % input: B, basis matrix; output: B_reduced, reduced basis matrix

        % Cong Ling, 2005

        % Based on the paper published later:

        % Ying Hung Gan, Cong Ling, and Wai Ho Mow, Complex lattice reduction 

        % algorithm for low-complexity full-diversity MIMO detection,

        % IEEE Trans. Signal Processing, vol. 57, pp. 2701-2710, July 2009. 



        M_ =  size(B_,2);         % number of columns; each column is a basis vector

        delta_ = .75;

        [Q_ R_] = qr(B_,0);          % use QR decomposition to do Gram-Schmit orthogonization

        beta = abs(diag(R_)).^2;      % squared length of orthogonal vectors (unnormalized)

        mu = R_ ./ (diag(diag(R_))*ones(M_));                          % normalize

        mu = mu.';              % to make it lower triangular

        k_ = 2;                  % k_ is the stage

        i_iteration = 0;

        while(i_iteration < 100*M_^2)                                    % loop until B_ cannot be reduced any more

            i_iteration = i_iteration + 1;                          % limit the maximum iteration number

            if (abs(real(mu(k_,k_-1))) > 0.5) | (abs(imag(mu(k_,k_-1))) > 0.5)   % to make it weakly reduced = abs(mu(k_,k_-1)) > 0.5

                [B_,mu] = size_reduce_k(B_,mu,k_,k_-1);

            end

            if(beta(k_) < (delta_ - abs(mu(k_,k_-1))^2) * beta(k_-1))    % swap the two if the k_-1th vector is in a sense longer than the kth

                b_ = B_(:,k_);

                B_(:,k_) = B_(:,k_-1);

                B_(:,k_-1) = b_;                                       % swap the kth column and the k_-1th column

                muswap = mu(k_-1,1:k_-2);

                mu(k_-1,1:k_-2) = mu(k_,1:k_-2);

                mu(k_,1:k_-2) = muswap;                               % update mu(k_-1,1:k_-2) and mu(k_,1:k_-2)

                old_muk = mu(k_+1:M_,k_);

                old_beta1 = beta(k_-1);

                old_betak = beta(k_);

                old_mu = mu(k_,k_-1);

                mu(k_+1:M_,k_) = mu(k_+1:M_,k_-1) - mu(k_+1:M_,k_) * mu(k_,k_-1);

                beta(k_-1) = beta(k_) + abs(mu(k_,k_-1))^2 * beta(k_-1); % update beta(k_-1),beta(k_), see the paper

                beta(k_) = beta(k_) * old_beta1 / beta(k_-1);

                mu(k_,k_-1) = mu(k_,k_-1)' * old_beta1 / beta(k_-1);

                mu(k_+1:M_,k_-1) = mu(k_+1:M_,k_-1) * mu(k_,k_-1) + old_muk * old_betak / beta(k_-1);

                if k_ > 2

                    k_ = k_-1;

                end

            else

                for i = k_-2 :-1: 1

                    if(abs(real(mu(k_,k_-1))) > 0.5) | (abs(imag(mu(k_,k_-1))) > 0.5)

                        [B_,mu] = size_reduce_k(B_,mu,k_,i);

                    end

                end

                if k_ < M_

                    k_ = k_ + 1;

                else

                    B_reduced = B_;

                    abs(mu);

                    return;

                end

            end

        end

        B_reduced = B_; % in i_iteration exceeds, simply returns B_ so far 'Warning: suboptimal CLLL basis'

        return;
        
    end



    function [B_,mu]=size_reduce_k(B_,mu,k_,j0)

        % in fact on step in size reduction

        eta = round(mu(k_,j0));

        B_(:,k_) = B_(:,k_) - eta * B_(:,j0);

        for i = 1 : j0-1

            mu(k_,i) = mu(k_,i) - eta * mu(j0,i);

        end

        mu(k_,j0) = mu(k_,j0) - eta;

        return;

    end


end