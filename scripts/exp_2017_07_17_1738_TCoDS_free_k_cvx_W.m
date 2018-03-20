% min_{k} |k-b| s.t. sum(k_i) = 4xi
% where b = V*sqrt(Dinv)*V^T*Ad*2/pi + 2/(pi|V|) * sum(Ad)
% and W = d0^Td0 = VDV^T is the (orthogonal) eigen decomposition of W
% see 19.7.17 in the notes.

%%%%%%%%%%%%%%
%% Settings %%
%%%%%%%%%%%%%%

%FNAME = 'sphere_s0.off';
%FNAME = 'sphere_s1.off';
%FNAME = 'ellipsoid_s4r.off';
%FNAME = 'round_cuber.off';
%FNAME = 'sphere_s0.off';
%FNAME = 'torus_fat_r2.off';
%FNAME = 'bunny.off';

%MESHES = {'sphere_s0.off', 'ellipsoid_s4r.off', 'round_cuber.off', 'torus_fat_r2.off', 'bunny.off', 'phands.off'};
MESHES = {'bunny.off'};

VERBOSE = true;
EPS = 1e-9;
EDGE_ALPHA = 0.2;
RESOLUTION = 1024;
AR = 1;

N_EIGENVALUES = 100;

theta0 = 0;

% T = datetime;
% OUT_FOLDER = fullfile('..', 'results', ...
%                       [num2str(T.Year), '.', sprintf('%02d', month(T)), '.', sprintf('%02d', T.Day)], ...
%                       [sprintf('%02d', T.Hour), sprintf('%02d', (T.Minute))]);
% if ~exist(OUT_FOLDER, 'dir')
%     mkdir(OUT_FOLDER);
% end

%%
for fname = MESHES
    disp(fname)
    % Find mesh file
    p = find_data_folder();
    fp = fullfile(p, fname{:});
    
    % Load mesh
    m = Mesh();
    m.loadTM(fp);
    scale = 2 * m.avg_length;
    
    [d0, d1] = get_exterior_derivatives(m);
    N = 4;          % Direction field degree (i.e., 4 for cross field)
    f0 = [1];       % Starting face
    %v0 = [0, 0, 1]; % Starting direction
    v0 = m.V(m.F(f0, 2), :) - m.V(m.F(f0, 1), :);
    
    Ad = get_gaussian_curvature(m);
    xi = 2 - 2*m.genus;
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% libigl MIQ (greedy rounding) %%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print_header('Solve with libigl MIQ')
    res1 = nrosy_wrapper(fp, f0, v0, N);
    r = res1.frame_diffs;
    local_frames = res1.local_frames;
%     
%     %%%%%%%%%%%%%%%%%%%%%%%
%     %% Solve with Gurobi %%
%     %%%%%%%%%%%%%%%%%%%%%%%
     print_header('Solve with Gurobi')
     W = d0'*d0;
     [V,D] = eig(full(W));
     
     [D, inds] = sort(diag(D));
     V = V(:, inds);
    
    Dinv = D;
    Dinv(Dinv < EPS) = inf;
    Dinv = 1 ./ Dinv;
    Dinv = diag(Dinv);
    D = diag(D);
%     
%     check_norm('V''*V', 'eye(size(V''*V))');
%     check_norm('V*V''', 'eye(size(V''*V))');
%     check_norm('D*Dinv', 'eye(size(D))');
%     check_norm('Dinv*D', 'eye(size(D))');
%     check_norm('V*D*V''', 'W');
%     
    Winv = V*Dinv*V';
%     
%     check_norm('W*Winv', 'eye(size(W))', 1-EPS);
%     
%     %%
%     b = V*sqrt(Dinv)*V'*(Ad)/(pi/2) + sum(Ad)/(pi/2)/m.nV;
%     %b = V*sqrt(Dinv)*V'*Ad/(pi/2);
%     br = round(b);
%     
%     cvx_solver gurobi
%     cvx_begin
%     cvx_solver_settings('TimeLimit', 10)
%         variable k_cvx(m.nV)
%         variable z(m.nV) integer
%         minimize ( norm(k_cvx - b) )
%         subject to
%             sum(k_cvx) == 4*xi
%             k_cvx == z
%     cvx_end
%     
%     %%
%     res2.k = k_cvx;
%     res2.u = Winv*(-Ad+(pi/2)*k_cvx);
%     res2.x = d0*res2.u;
%     [res2.nrosy, res2.theta] = create_direction_field(m, res2.x, f0, theta0, N);
%     res2.S_inds = find(abs(k_cvx) > EPS);
%     res2.E = norm(res2.x)^2;
%     
%     check_norm('d0''*res2.x', '-Ad+(pi/2)*k_cvx');
%     
%     inds = find(abs(k_cvx) > EPS);
%     I = [inds, k_cvx(inds)/N];
%     [res3.x, res3.E] = create_trivial_connection(m, I, N, VERBOSE);
%     [res3.nrosy, res3.theta] = create_direction_field(m, res3.x, f0, theta0, N);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Solve with closest lattice %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %b = [V*sqrt(Dinv)*V'*(Ad)/(pi/2) + sum(Ad)/(pi/2)/m.nV; -4*xi-1];
    b = [V*sqrt(Dinv)*V'*(Ad)/(pi/2); -4*xi];
    %b = [sqrt(D)*Dinv*V'*Ad*2/pi; -4*xi];
    b = b - sum(b)/(length(b));
    check_norm('sum(b)', '0');
    br = round(b);
    
    %b = [0.168, 0.462, -0.631];
    %br = round(b);
    
    while 1
        %br = round(b);   
        
        def = sum(br);
        delta = b - br; %br - b;
        [delta, inds] = sort(delta);
        b = b(inds);
        br = br(inds);
        
        if def == 0
            br(inds) = br;
            b(inds) = b;
            break;
        elseif def > 0
            br(1:def) = br(1:def-1) - 1;
        else
            br(end+def+1:end) = br(end+def+1:end) + 1;
        end
        
        br(inds) = br;
        b(inds) = b;
        
    end
    
    %res2.z = V*sqrt(D)*V'*br(1:end-1);
    %res2.k = br(1:end-1);
    tic
    res2.k = round(V*sqrt(D)*V'*br(1:end-1));
    
    for i = 1:m.nV
    if abs(res2.k) > 0
        continue
    end
    disp(i)
    res2.k(i) = 1;
    
    res2.u = Winv*(-Ad+(pi/2)*res2.k);
    res2.x = d0*res2.u;
    %[res2.nrosy, res2.theta] = create_direction_field(m, res2.x, f0, theta0, N);
    %res2.S_inds = find(abs(res2.k) > EPS);
    res2.E = norm(res2.x)^2;
    E(i) = res2.E;
    
    res2.k(i) = 0;
    
    end
    toc
    
    check_norm('d0''*res2.x', '-Ad+(pi/2)*res2.k');
    
    inds = find(abs(res2.k) > EPS);
    I = [inds, res2.k(inds)/N];
    [res3.x, res3.E] = create_trivial_connection(m, I, N, VERBOSE);
    [res3.nrosy, res3.theta] = create_direction_field(m, res3.x, f0, theta0, N);
    
    [res4.theta, res4.p] = TCoDS_to_MIQ(res2.x, res2.k, r, theta0, Ad, d0, d1, res1.period_is_fixed);
    res4.E = E_MIQ(m, res4.theta, r, res4.p, N);
    
    [res5.x, res5.k] = MIQ_to_TCoDS(res1.theta, res1.p, r, Ad, d0, d1);
    res5.nrosy = create_direction_field(m, res5.x, f0, theta0, N, local_frames, r);
    res5.E = norm(res5.x)^2;
    
    
    %res2.E1 = norm(V*sqrt(Dinv)*V'*(-Ad+(pi/2)*res2.k))^2;
    %res2.E2 = norm((pi/2)*res2.k-V*sqrt(Dinv)*V'*Ad)^2;
    %res2.E3 = norm(res2.k-(2/pi)*V*sqrt(Dinv)*V'*Ad)^2;
    %res2.E4 = norm(res2.k-(2/pi)*V*sqrt(Dinv)*V'*Ad)^2 + norm(br(end)-(-4*xi))^2;
    
    res2.E1 = norm(res2.k - (2/pi)*V*sqrt(Dinv)*V'*Ad)^2;
    res2.E2 = norm(V*sqrt(Dinv)*V'*(-Ad+(pi/2)*res2.k))^2;
    
    res5.E1 = norm(res5.k - (2/pi)*V*sqrt(Dinv)*V'*Ad)^2;
    res5.E2 = norm(V*sqrt(Dinv)*V'*(-Ad+(pi/2)*res5.k))^2;
    
    
    
    %%
    cvx_solver gurobi
    cvx_begin
    %cvx_solver_settings('TimeLimit', 10)
        variable x(m.nV)
        variable z(m.nV) integer
        minimize( norm(x - (2/pi)*V*sqrt(Dinv)*V'*Ad ) + norm(x) )
        subject to
            sum( V*sqrt(D)*V'*x ) == 4*xi
            %x == z 
    cvx_end
    
    
    
    %%
%     [U,S,V] = svds(W,m.nV);
%     Sinv = diag(S);
%     Sinv(Sinv < EPS) = inf;
%     Sinv = 1 ./ Sinv;
%     Sinv = diag(Sinv);
%     Winv = V * Sinv * V';
%     
%     check_norm('U*S*V''', 'W');
%     check_norm('U*U''', 'eye(m.nV)');
%     check_norm('U''*U', 'eye(m.nV)');
%     check_norm('V''*V', 'eye(m.nV)');
%     check_norm('U', 'V');
%     check_norm('S*Sinv', 'eye(size(S))');
%     %check_norm('S*inv(S)', 'eye(size(S))');
%     check_norm('W*Winv', 'eye(size(W))');
%     
%     br = round(V*sqrt(Sinv)*V'*(-Ad)/(pi/2));
%     
%     %%
%     cvx_solver gurobi
%     cvx_begin
%     cvx_solver_settings('TimeLimit', 10)
%         variable k_cvx(m.nV)
%         variable z(m.nV) integer
%         minimize ( norm(k_cvx-(res1.S~=0)) )
%         subject to
%             sum(k_cvx) == 4*xi
%             k_cvx == z
%     cvx_end
%     
%     u = Winv*(-Ad+(pi/2)*k_cvx);
%     x = d0*u;
%     
%     check_norm('d0''*x', '-Ad+(pi/2)*k_cvx');
%     
%     inds = find(abs(k_cvx) > EPS);
%     I = [inds, k_cvx(inds)/N];
%     [res3.x, res3.E] = create_trivial_connection(m, I, N, VERBOSE);
    
%     tic
%     %P = [d0'*d0, zeros(m.nV); zeros(m.nV), eye(m.nV)];
%     P = [d0'*d0, zeros(m.nV); zeros(m.nV), zeros(m.nV)];
%     q = zeros(2*m.nV, 1);
%     r = 0;
%     A = [d0'*d0, -(pi/2)*eye(m.nV);
%          zeros(1, m.nV), ones(1, m.nV);
%          ones(1, m.nV), zeros(1, m.nV)];
%     b = [-Ad; 4*xi; 0];
%     cvx_solver gurobi
%     cvx_begin
%     cvx_solver_settings('TimeLimit', 10)
%         variable x_cvx(2*m.nV)
%         variable z(m.nV) binary
%         minimize (0.5*x_cvx'*P*x_cvx)
%         subject to
%             A*x_cvx == b
%             x_cvx(m.nV+1:2*m.nV) == z
%             x_cvx(m.nV+1:2*m.nV) <= 10
%             x_cvx(m.nV+1:2*m.nV) >= -10
%             %sum(x_cvx(m.nV+1:2*m.nV)) == N*xi
%             %sum(x_cvx(1:m.nV)) == 0
%             %ones(m.nV, 1)'*x_cvx(m.nV+1:2*m.nV) <= 10
%             %ones(m.nV, 1)'*x_cvx(m.nV+1:2*m.nV) >= -10
%     cvx_end
%     toc
%     
%     res2.u = x_cvx(1:m.nV);
%     res2.x = d0*res2.u;
%     res2.k = x_cvx(m.nV+1:end);
%     [res2.nrosy, res2.theta] = create_direction_field(m, res2.x, f0, theta0, N);
%     res2.E = norm(res2.x)^2;
%     res2.S_inds = find(res2.k > EPS);
%     
%     check_norm('d0''*res2.x-(pi/2)*res2.k','-Ad');
%     
%     res2.E
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %% Run TCoDS with k from previous step %%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     inds = res2.S_inds;
%     I = [inds, res2.k(inds)/N];
%     %theta0 = [res1.local_frames(1,:); res1.local_frames(m.nF+1,:)] * v0';
%     %theta0 = atan2(theta0(2), theta0(1));
%     %theta0 = 0;
%     [res3.x, res3.E] = create_trivial_connection(m, I, N, VERBOSE);
%     [res3.nrosy, res3.theta] = create_direction_field(m, res3.x, f0, theta0, N, res1.local_frames, res1.frame_diffs); 
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%
%     %% CoV from cvx to MIQ %%
%     %%%%%%%%%%%%%%%%%%%%%%%%%
%     [res4.theta, res4.p] = TCoDS_to_MIQ(res2.x, res2.k, res1.frame_diffs, theta0, Ad, d0, d1, res1.period_is_fixed);
%     res4.E = E_MIQ(m, res4.theta, res1.frame_diffs, res4.p, N);
%     res4.nrosy = angles_to_nrosy(res4.theta, res1.local_frames, N);
%     
%     %%%%%%%%%%%%%%
%     %% MIQ ADMM %%
%     %%%%%%%%%%%%%%
%     %      - k - | -- u --  
%     % P = [  0   | d0'*d0 ]
%     %     [  0   | 0      ]
%     %
%     %     [ -pi/2*I | d0'*d0 ]
%     % A = [ 0       | 0      ]
%     %     [ 1       | 0      ]
%     %
%     %     [-Ad   ]
%     % b = [ 4*xi ]
%     %     [ 0    ]
%     xtmp = x_cvx;
%     xtmp(1:m.nV) = x_cvx(m.nV+1:end);
%     xtmp(m.nV+1:end) = x_cvx(1:m.nV);
%     P_admm = ([zeros(m.nV), d0'*d0; 
%                zeros(m.nV), zeros(m.nV)]);
%     q_admm = zeros(2*m.nV, 1);
%     A_admm = [-(pi/2)*eye(m.nV), d0'*d0;
%               ones(1, m.nV), zeros(1, m.nV);
%               zeros(1, m.nV), ones(1, m.nV)];
%     r_admm = 0;
%     b_admm = [-Ad; 4*xi; 0];
%     Problem.objective.P = full(P_admm);
%     Problem.objective.q = q_admm;
%     Problem.objective.r = r_admm;
%     %Problem.constraint.A = full(A);
%     %Problem.constraint.b = b;
%     Problem.constraint.A = full(A_admm);
%     Problem.constraint.b = b_admm;
%     Problem.constraint.l1 = m.nV;
%     Problem.constraint.l2 = 0;
%     Problem.constraint.l3 = 0;
%     Problem.constraint.l4 = 0;
%     
%     maxiter = 5000; % max number of iterations
%     repeat = 10; % number of initializations
%     xs = solver_miqp_admm(Problem, 1, maxiter, repeat);
% 
%     f_admm = Inf;
%     x_admm = zeros(2*m.nV,1);
%     res_thr = 1e-3; % feasibilty threshold 
%     thr = zeros(repeat*maxiter, 1);
%     f = zeros(repeat*maxiter, 1);
%     i = 1;
%     normA = norm(full(A_admm));
%     for j=1:repeat
%         for k=1:maxiter
%             x = xs(:,(j-1)*maxiter+k);
%             thr(i) = norm(A_admm*x-b_admm);
%             f(i) = 0.5*x'*P_admm*x + q_admm'*x + r_admm;
%             %fprintf('%f, %f\n', thr(i), f(i))
%             if thr(i) < res_thr
%                 %f = 0.5*x'*P*x + q'*x + r;
%                 if f(i) < f_admm
%                     f_admm = f(i);
%                     x_admm = x;
%                 end
%             end
%             i = i + 1;
%         end
%     end
%     
%     %assert(f_admm < Inf);
%     min(thr)
%     fprintf('The best value found by ADMM is %3s\n', f_admm);
%     
%     res5.u = x_admm(m.nV+1:end);
%     res5.k = x_admm(1:m.nV);
%     res5.x = d0*res5.u;
%     
%     check_norm('d0''*res5.x-(pi/2)*res5.k','-Ad')
% 
%     %assert(norm(res4.x - res2.x) < EPS)
%     %norm(d1*res4.x)
    
    %%%%%%%%%%%
    %% Plots %%
    %%%%%%%%%%%
    nPlots = 2;
    mPlots = 2;
    spi = 1;
    figi = 1;
    ha = zeros(nPlots*mPlots,1);
    scale = 2*m.avg_length;

    figure(figi); clf(figi)
    ha(1) = subplot(nPlots, mPlots, spi);
    spi = spi + 1;
    hold on
    draw_constraints(m, f0, v0)
    draw_direction_field(m, res1.nrosy, res1.S_inds, scale)
    hold off
    title({['libigl MIQ (greedy rounding)'], ['E = ', num2str(res1.E)]})

    ha(2) = subplot(nPlots, mPlots, spi);
    spi = spi + 1;
    hold on
    draw_constraints(m, f0, v0);
    draw_direction_field(m, res2.nrosy, res2.S_inds, scale)
    hold off
    title({'TCoDS free k (CVX)', ['E = ', num2str(res2.E)]})

    ha(3) = subplot(nPlots, mPlots, spi);
    spi = spi + 1;
    hold on
    draw_constraints(m, f0, v0);
    draw_direction_field(m, res3.nrosy, res2.S_inds, scale)
    hold off
    title({'TCoDS with CVX singularities', ['E = ', num2str(res3.E)]})
    
    ha(4) = subplot(nPlots, mPlots, spi);
    spi = spi + 1;
    m.draw(b(1:end-1))
    title('V*sqrt(Dinv)*V''*(Ad)/(pi/2)');
    colorbar
    hold on
    % plot singularities
    I = res2.S_inds;
    for j = 1:size(I,1)
        fid = I(j, 1);
        H = plot3( m.V(fid,1), m.V(fid,2), m.V(fid,3), 'r.' );
        set( H, 'MarkerSize', 40 );
    end
    hold off
    %subplot(212)
    %plot(diag(D), '.')
    %title('diag(D)')
    
    

%     ha(4) = subplot(nPlots, mPlots, spi);
%     spi = spi + 1;
%     hold on;
%     draw_constraints(m, f0, v0);
%     draw_direction_field(m, res4.nrosy, res2.S_inds, scale)
%     hold off
%     title({'TCoDScvx(x,k) --> MIQ(p,theta)', ['E = ', num2str(res4.E)]})

    hlink = linkprop(ha, {'CameraPosition', 'CameraUpVector'});
    figi = figi + 1;
end