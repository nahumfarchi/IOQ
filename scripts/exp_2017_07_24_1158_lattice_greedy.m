function [] = exp_2017_07_24_1158_lattice_greedy()
ME = [];
try
    
    ABOUT = ['minimize norm(M*k-t)^2 \r\n', ...
             's.t. sum(k_i) = m \r\n', ...
             '\r\n', ...
             'where M = sqrt(Dinv)*V \r\n', ...
             '\tb = (2/pi)*Ad \r\n', ...
             '\tt = (2/pi)*M*Ad \r\n', ...
             '\tm = 4*xi \r\n', ...
             'using greedy lattice with {+1,-1} pair search\r\n\r\n'];

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
    SAVE = false;

    if PLOT && SAVE
        OUT_FOLDER = create_time_stamped_folder(fullfile('..', 'results'));
        LOG = fopen(fullfile(OUT_FOLDER, 'log.txt'), 'w')
    else
        LOG = -1;
    end

    log_and_print(LOG, ABOUT);
    log_and_print(LOG, 'logid : %d\r\n', LOG);
    %%
    for fname = MESHES
        log_and_print(LOG, 'Loading %s...\r\n', fname{:});
        p = find_data_folder();
        fp = fullfile(p, fname{:});

        mesh = Mesh();
        mesh.loadTM(fp);
        scale = 2 * mesh.avg_length;

        Ad = get_gaussian_curvature(mesh);
        xi = 2 - 2*mesh.genus;
        [d0, d1] = get_exterior_derivatives(mesh);

        % TODO doesn't work with N!=0 for some reason
        N = 4;    % Direction field degree (i.e., 4 for cross field)
        f0 = [1]; % Starting face
        v0 = mesh.V(mesh.F(f0, 2), :) - mesh.V(mesh.F(f0, 1), :);

        log_and_print(LOG, 'N : %d\r\n', N);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % libigl MIQ (greedy rounding) %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if RUN_LIBIGL_MIQ
            log_and_print(LOG, '\r\nSolving with libigl MIQ...\r\n');
            [T, res1] = evalc('nrosy_wrapper(fp, f0, v0, N)');
            log_and_print(LOG, T);
        end

        %%%%%%%%%%%%%%%%%%
        % Greedy lattice %
        %%%%%%%%%%%%%%%%%%
        %%
        log_and_print(LOG, '\r\nSolving with greedy lattice...\r\n');

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

        M = sqrt(Dinv)*V';
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
        check_norm('sum(k)', '4*xi', 'Log', LOG);
        Mcols = sum(M, 2); % Keep track of sum(M_i) over S
        tic
        Mnorms = norms(M, [], 1).^2;
        krnd = zeros(size(k));
        kfix = zeros(size(k));
        MM = M'*M;
        tM = t'*M;
        VV = mesh.VVAdj;
        free_vars = (1:mesh.nV)';
        log_and_print(LOG, 'Rounding variables...\r\n');
        while nfree > 1
            perc = 100*(length(k)-nfree)/length(k);
            nTenPerc = floor(length(k) / 10);
            if mod(length(k)-nfree, nTenPerc) == 0
            %if floor(mod(perc, 10)) == 0
            %if mod(nfree, 100) == 0
                %disp([num2str(100*(length(k)-nfree)/length(k)), '%...'])
                %log_and_print(LOG, '%g%%\r\n', 100*(length(k)-nfree)/length(k));
                log_and_print(LOG, '%g%%...\r\n', perc);
            end

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
                [best_rounding_E, best_rounding_z, best_rounding_q] = ...
                    find_best_z(M(:,i), Mcols, Mkt, k(i), m, xi, nfree, integer_candidates);
                Mcols = Mcols + M(:, i);

                if best_rounding_E < best_E
                    best_i = i;
                    best_E = best_rounding_E;
                    best_z = best_rounding_z;
                    best_q = best_rounding_q;
    %                 best_inds = best_rounding_inds;
                end           
            end

            Mcols = Mcols - M(:, best_i);
            S(best_i) = 0;
            free_vars = setdiff(free_vars, best_i);
            nfree = nfree - 1;
            k(best_i) = best_z;
            k(S) = k(S) + best_q;
        end

        k(end) = m - sum(k(1:end-1));
        k = round(k);
        check_norm('sum(k)', '4*xi', 'Log', LOG);

        log_and_print(LOG, 'Searching for pairs to set to +1, -1...\r\n');
        C = nchoosek(1:length(k), 2);
        Mk = M*k;
        for i = 1:size(C,1)
            perc = 100 * i / size(C,1);
            nTenPerc = floor(size(C,1) / 10);
            if mod(i, nTenPerc) == 0
            %if floor(mod(perc, 10)) == 0
            %if mod(i, 1000) == 0
                %log_and_print(LOG, '%s%%...\r\n', num2str(100 * i / size(C,1)));
                log_and_print(LOG, '%g%%...\r\n', perc);
            end
            E_old = norm(Mk-t)^2;

            ktmp1 = k;
            j1 = C(i, 1);
            j2 = C(i, 2);
            ktmp1(j1) = ktmp1(j1) + 1;
            ktmp1(j2) = ktmp1(j2) - 1;
            E1 = norm(Mk + M(:,j1) - M(:,j2) - t)^2;

            ktmp2 = k;
            j1 = C(i, 1);
            j2 = C(i, 2);
            ktmp2(j1) = ktmp2(j1) - 1;
            ktmp2(j2) = ktmp2(j2) + 1;
            E2 = norm(Mk - M(:,j1) + M(:,j2) - t)^2;

            best_E = min(min(E1, E2), E_old);
            if abs(best_E - E1) < EPS
                k = ktmp1;
                Mk = Mk + M(:,j1) - M(:,j2);
            elseif abs(best_E - E2) < EPS
                k = ktmp2;
                Mk = Mk - M(:,j1) + M(:,j2);
            end

        end

        check_norm('sum(k)', '4*xi', 'Log', LOG);

        lattice_time = toc;
        log_and_print(LOG, 'Greedy lattice time: %g\r\n', lattice_time);

        res2.k = k;
        res2.u = Winv*(-Ad+(pi/2)*res2.k);
        res2.x = d0*res2.u;
        res2.E = norm(res2.x)^2;
        %res2.nrosy = TCODS(mesh, S, degree, f0, theta0, VERBOSE);
        %[res2.nrosy, res2.theta] = create_direction_field(mesh, res2.x, f0, theta0, N);
        [res2.ffield, res2.theta]  = connection_to_nrosy(mesh, res2.x, f0, theta0, N);
        res2.S_inds = find(abs(res2.k) > EPS);
        res2.S = [res2.S_inds, res2.k(res2.S_inds)/N];

        check_norm('d0''*res2.x', '-Ad+(pi/2)*res2.k', 'Log', LOG);

        %inds = find(abs(res2.k) > EPS);
        %I = [inds, res2.k(inds)/N];
        res3 = TCODS(mesh, res2.S, f0, theta0, N, VERBOSE);
        %[res3.x, res3.E] = create_trivial_connection(mesh, res2.S, N, VERBOSE);
        %[res3.ffield, res3.theta] = connection_to_nrosy(mesh, res3.x, f0, theta0, N);

        [res4.theta, res4.p] = TCoDS_to_MIQ(res3.x, res2.k, res1.frame_diffs, theta0, Ad, d0, d1, res1.period_is_fixed);
        res4.E = E_MIQ(mesh, res4.theta, res1.frame_diffs, res4.p, N);
        res4.ffield = angles_to_nrosy(res4.theta, res1.local_frames, N);

        %%%%%%%%%
        % Plots %
        %%%%%%%%%
        %%
        if ~PLOT
            continue
        end

        log_and_print(LOG, '\r\nPlotting...\r\n');

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
        draw_direction_field(mesh, res1.ffield, res1.S_inds, scale)
        hold off
        title({['libigl MIQ (greedy rounding)'], ['E = ', num2str(res1.E)]})

        ha(2) = subplot(nPlots, mPlots, spi);
        spi = spi + 1;
        hold on
        draw_constraints(mesh, f0, v0);
        draw_direction_field(mesh, res2.ffield, res2.S_inds, scale)
        hold off
        title({'Greedy lattice', ['E = ', num2str(res2.E)]})

        ha(3) = subplot(nPlots, mPlots, spi);
        spi = spi + 1;
        hold on
        draw_constraints(mesh, f0, v0);
        draw_direction_field(mesh, res3.ffield, res2.S_inds, scale)
        hold off
        title({'TCoDS with greedy lattice singularities', ['E = ', num2str(res3.E)]})

        ha(4) = subplot(nPlots, mPlots, spi);
        spi = spi + 1;
        hold on;
        draw_constraints(mesh, f0, v0);
        draw_direction_field(mesh, res4.ffield, res2.S_inds, scale)
        hold off
        title({'TCoDS(x,k) --> MIQ(p,theta)', ['E = ', num2str(res4.E)]})

        hlink = linkprop(ha, {'CameraPosition', 'CameraUpVector'});
        figi = figi + 1;

        if SAVE
            log_and_print(LOG, '\r\nSaving plots...\r\n');

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

        log_and_print(LOG, '\r\nDone.\r\n');
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


%%%%%%%%%%%%%%%%%%%%
% Helper Functions %
%%%%%%%%%%%%%%%%%%%%

    function [Eopt, zopt, qopt] = find_best_z(Mi, Mcols, Mkt, ki, m, xi, nfree, integer_candidates)
        % Find the best integer to round to
        Eopt = inf;
        for z = integer_candidates
            Mkrnd = Mi*(z - ki);
            q = (m - (z-ki+4*xi)) / (nfree - 1);
            Mkfix = q*Mcols;

            % 0.635s, E=3.7287 on sphere_s0 
            rounding_E = norm(Mkrnd+Mkfix)^2 + 2*Mkt'*(Mkrnd + Mkfix);

%                 inds = VV(i,:)';
%                 inds(inds~=0 & S==0) = 0;
%                 inds(i) = 0;
%                 q = (m - (z-k(i)+4*xi)) / sum(inds);
%                 
%                 
%                 
% %                 %rounding_E = norm(Mkrnd)^2+norm(Mkfix)^2+2*Mkrnd'*Mkfix + 2*Mkt'*(Mkrnd + Mkfix);
% %                 %rounding_E = norm(Mkrnd)^2+norm(Mkfix)^2+2*Mkrnd'*Mkfix + ...
% %                 %    2*k'*M'*(Mkrnd+Mkfix) - 2*t'*(Mkrnd+Mkfix);
% %                 %rounding_E = (z-k(i))^2 * Mnorms(i);
% %                 krnd(i) = z-k(i);
% %                 S(i) = 0;
% %                 %inds = VV(i, :);
% %                 %inds = randsample(setdiff(find(S~=0), i), length(setdiff(find(S~=0), i)));
% %                 %inds = datasample(find(S~=0), min(10, nfree-1), 'Replace', false);
% %                 %inds = 1:10;
% %                 inds = VV(i,:);
% %                 %inds = datasample(free_vars, min(10, nfree-1), 'Replace', false);
% %                 %kfix(inds) = q;
% %                 q = (m - (z-k(i)+4*xi)) / sum(inds);
% %                 
% %                 %kfix(S) = q;
% %                 %rounding_E = rounding_E + q*norm(Mcols)^2 + 2*krnd'*MM*kfix;
% %                 %rounding_E = rounding_E + 2*Mkt'*(Mkrnd + Mkfix);
%                 zki = z - k(i);
%                 rounding_E = zki^2*Mnorms(i);
%                 rounding_E = rounding_E + q^2*norm(M(:,inds))^2;
%                 rounding_E = rounding_E + 2*zki*sum(MM(i,inds))*q;
%                 rounding_E = rounding_E + 2*kMM(i)*zki;
%                 rounding_E = rounding_E + 2*sum(kMM(inds))*q;
%                 rounding_E = rounding_E - 2*tM(i)*zki;
%                 rounding_E = rounding_E - 2*sum(tM(inds))*q;
%                 
%                 ktmp = k;
%                 ktmp(i) = z;
%                 ktmp(inds) = ktmp(inds) + q;
%                 assert(check_norm('sum(ktmp)', '4*xi'));
% %                     %2*tM*(krnd+kfix);
% %                     %2*kMM*(krnd+kfix) - ...
% %                     %2*tM*(krnd+kfix);
% %                     %2*Mkt'*(Mkrnd + Mkfix);
% %                 S(i) = 1;
% %                 %kfix(inds) = 0;
% %                 krnd(i) = 0;              

            if rounding_E < Eopt
                Eopt = rounding_E;
                zopt = z;
                qopt = q;
%                     best_rounding_inds = inds;
            end
        end
    end

end