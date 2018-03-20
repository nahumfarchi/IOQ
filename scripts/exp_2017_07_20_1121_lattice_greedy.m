% minimize norm(M*k-t)^2
% s.t. sum(k_i) = m
%
% where M = sqrt(Dinv)*V'
%       t = (2/pi)*M*Ad
%       m = 4*xi

%%%%%%%%%%%%
% Settings %
%%%%%%%%%%%%
%%
%MESHES = {'sphere_s0.off', 'ellipsoid_s4r.off', 'round_cuber.off', 'torus_fat_r2.off', 'bunny.off', 'phands.off'};
%MESHES = {'sphere_s0.off', 'bumpy.off', 'round_cuber.off', 'cow.off', 'bunny.off', 'phands.off'};
MESHES = {'sphere_s0.off'};

VERBOSE = true;
EPS = 1e-9;
EDGE_ALPHA = 0.2;
RESOLUTION = 1024;
AR = 1;

theta0 = 0;

RUN_LIBIGL_MIQ = true;
PLOT = false;

if PLOT
    OUT_FOLDER = create_time_stamped_folder(fullfile('..', 'results'));
end

%%
for fname = MESHES
    disp(fname)
    p = find_data_folder();
    fp = fullfile(p, fname{:});
    
    mesh = Mesh();
    mesh.loadTM(fp);
    scale = 2 * mesh.avg_length;

    Ad = get_gaussian_curvature(mesh);
    xi = 2 - 2*mesh.genus;
    [d0, d1] = get_exterior_derivatives(mesh);
    
    N = 1;    % Direction field degree (i.e., 4 for cross field)
    f0 = [1]; % Starting face
    v0 = mesh.V(mesh.F(f0, 2), :) - mesh.V(mesh.F(f0, 1), :);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % libigl MIQ (greedy rounding) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if RUN_LIBIGL_MIQ
        print_header('Solve with libigl MIQ')
        res1 = nrosy_wrapper(fp, f0, v0, N);
    end
    
    %%%%%%%%%%%%%%%%%%
    % Greedy lattice %
    %%%%%%%%%%%%%%%%%%
    %%
    print_header('Solve with greedy lattice')
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
    
    M = V*sqrt(Dinv)*V';
    b = (2/pi)*Ad;
    t = (2/pi)*M*Ad;
    %M = V*D*V';
    %t = V*D*sqrt(D)*t;
    m = 4*xi;
    
    k = b;
    S = logical(ones(size(k))); % Indices of non-integer variables
    best_rounding = nan(size(k));
    nfree = mesh.nV;
    integer_candidates = [2, 1, 0, -1, -2];
    check_norm('sum(k)', '4*xi');
    Mcols = sum(M, 2); % Keep track of sum(M_i) over S
    tic
    Mnorms = norms(M, [], 1).^2;
    krnd = zeros(size(k));
    kfix = zeros(size(k));
    MM = M'*M;
    tM = t'*M;
    VV = mesh.VVAdj;
    free_vars = (1:mesh.nV)';
    while nfree > 1
        disp(nfree)
        best_i = 1;
        best_E = inf;
        Mkt = M*k - t;
        kMM = k'*MM;
        
        % Find the best ki to round
        for i = 1:length(S)
            if ~S(i)
                continue
            end
            
            Mcols = Mcols - M(:, i);
            best_rounding_E = inf;
            
            % Find the best integer to round to
            for z = integer_candidates
                Mkrnd = M(:,i)*(z - k(i));
                q = (m - (z-k(i)+4*xi)) / (nfree - 1);
                Mkfix = q*Mcols;
                
                % 0.635s, E=3.7287 on sphere_s0 
                rounding_E = norm(Mkrnd+Mkfix)^2 + 2*Mkt'*(Mkrnd + Mkfix);
                
                
                
                
%                 %rounding_E = norm(Mkrnd)^2+norm(Mkfix)^2+2*Mkrnd'*Mkfix + 2*Mkt'*(Mkrnd + Mkfix);
%                 %rounding_E = norm(Mkrnd)^2+norm(Mkfix)^2+2*Mkrnd'*Mkfix + ...
%                 %    2*k'*M'*(Mkrnd+Mkfix) - 2*t'*(Mkrnd+Mkfix);
%                 %rounding_E = (z-k(i))^2 * Mnorms(i);
%                 krnd(i) = z-k(i);
%                 S(i) = 0;
%                 %inds = VV(i, :);
%                 %inds = randsample(setdiff(find(S~=0), i), length(setdiff(find(S~=0), i)));
%                 %inds = datasample(find(S~=0), min(10, nfree-1), 'Replace', false);
%                 %inds = 1:10;
%                 inds = VV(i,:);
%                 %inds = datasample(free_vars, min(10, nfree-1), 'Replace', false);
%                 %kfix(inds) = q;
%                 q = (m - (z-k(i)+4*xi)) / sum(inds);
%                 
%                 %kfix(S) = q;
%                 %rounding_E = rounding_E + q*norm(Mcols)^2 + 2*krnd'*MM*kfix;
%                 %rounding_E = rounding_E + 2*Mkt'*(Mkrnd + Mkfix);
%                 zki = z - k(i);
%                 rounding_E = zki^2*Mnorms(i);
%                 rounding_E = rounding_E + q^2*norm(Mcols)^2;
%                 rounding_E = rounding_E + 2*zki*sum(MM(i,inds))*q;
%                 rounding_E = rounding_E + 2*kMM(i)*zki;
%                 rounding_E = rounding_E + 2*sum(kMM(inds))*q;
%                 rounding_E = rounding_E - 2*tM(i)*zki;
%                 rounding_E = rounding_E - 2*sum(tM(inds))*q;
%                     %2*tM*(krnd+kfix);
%                     %2*kMM*(krnd+kfix) - ...
%                     %2*tM*(krnd+kfix);
%                     %2*Mkt'*(Mkrnd + Mkfix);
%                 S(i) = 1;
%                 %kfix(inds) = 0;
%                 krnd(i) = 0;              
                
                if rounding_E < best_rounding_E
                    best_rounding_E = rounding_E;
                    best_rounding_z = z;
                    best_rounding_q = q;
                end
            end
            
            if best_rounding_E < best_E
                best_i = i;
                best_E = best_rounding_E;
                best_z = best_rounding_z;
                best_q = best_rounding_q;
                %best_k = best_rounding;
            end
            
            Mcols = Mcols + M(:, i);
        end
        
        Mcols = Mcols - M(:, best_i);
        S(best_i) = 0;
        free_vars = setdiff(free_vars, best_i);
        nfree = nfree - 1;
        k(best_i) = best_z;
        k(S) = k(S) + best_q;
    end
    toc
    
    k(end) = m - sum(k(1:end-1));
    k = round(k);
    check_norm('sum(k)', '4*xi');
    
    res2.k = k;
    res2.u = Winv*(-Ad+(pi/2)*res2.k);
    res2.x = d0*res2.u;
    res2.E = norm(res2.x)^2;
    [res2.nrosy, res2.theta] = create_direction_field(mesh, res2.x, f0, theta0, N);
    res2.S_inds = find(abs(res2.k) > EPS);
    
    check_norm('d0''*res2.x', '-Ad+(pi/2)*res2.k');
    
    inds = find(abs(res2.k) > EPS);
    I = [inds, res2.k(inds)/N];
    [res3.x, res3.E] = create_trivial_connection(mesh, I, N, VERBOSE);
    [res3.nrosy, res3.theta] = create_direction_field(mesh, res3.x, f0, theta0, N);
    
    [res4.theta, res4.p] = TCoDS_to_MIQ(res3.x, res2.k, res1.frame_diffs, theta0, Ad, d0, d1, res1.period_is_fixed);
    res4.E = E_MIQ(mesh, res4.theta, res1.frame_diffs, res4.p, N);
    res4.nrosy = angles_to_nrosy(res4.theta, res1.local_frames, N);
    
    %%%%%%%%%
    % Plots %
    %%%%%%%%%
    %%
    if ~PLOT
        continue
    end
    
    nPlots = 2;
    mPlots = 2;
    spi = 1;
    figi = 1;
    ha = zeros(nPlots*mPlots,1);
    scale = 2*mesh.avg_length;

    figure(figi); clf(figi)
    ha(1) = subplot(nPlots, mPlots, spi);
    spi = spi + 1;
    hold on
    draw_constraints(mesh, f0, v0)
    draw_direction_field(mesh, res1.nrosy, res1.S_inds, scale)
    hold off
    title({['libigl MIQ (greedy rounding)'], ['E = ', num2str(res1.E)]})

    ha(2) = subplot(nPlots, mPlots, spi);
    spi = spi + 1;
    hold on
    draw_constraints(mesh, f0, v0);
    draw_direction_field(mesh, res2.nrosy, res2.S_inds, scale)
    hold off
    title({'Greedy lattice', ['E = ', num2str(res2.E)]})

    ha(3) = subplot(nPlots, mPlots, spi);
    spi = spi + 1;
    hold on
    draw_constraints(mesh, f0, v0);
    draw_direction_field(mesh, res3.nrosy, res2.S_inds, scale)
    hold off
    title({'TCoDS with greedy lattice singularities', ['E = ', num2str(res3.E)]})
    
    ha(4) = subplot(nPlots, mPlots, spi);
    spi = spi + 1;
    hold on;
    draw_constraints(mesh, f0, v0);
    draw_direction_field(mesh, res4.nrosy, res2.S_inds, scale)
    hold off
    title({'TCoDS(x,k) --> MIQ(p,theta)', ['E = ', num2str(res4.E)]})

    hlink = linkprop(ha, {'CameraPosition', 'CameraUpVector'});
    figi = figi + 1;
    
    [~, name_only, ~] = fileparts(fname{:});
    pngname = fullfile(OUT_FOLDER, [name_only, '.png']);
    %pngname = sprintf('%s_%04d.png',filename,i);
                
    dpi = get(0, 'ScreenPixelsPerInch');
    in = RESOLUTION/dpi;

    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 AR*in in];
    fig.PaperPositionMode = 'manual';
    print(pngname, '-dpng', '-r0')
end