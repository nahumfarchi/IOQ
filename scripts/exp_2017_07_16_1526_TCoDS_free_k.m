% Solve argmin_{u,k} \|d_0u\|_2^2
%       s.t. d_0^Td_0u = -A_d + (pi/2)k
%            k \in \mathbb{Z}
%            \sum k(i) = \sum A_d(i) / (pi/2)
%            \sum u(i) = 0

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
MESHES = {'bumpy.off'};

VERBOSE = true;
EPS = 1e-9;
EDGE_ALPHA = 0.2;
RESOLUTION = 1024;
AR = 1;

theta0 = 0;

% T = datetime;
% OUT_FOLDER = fullfile('..', 'results', ...
%                       [num2str(T.Year), '.', sprintf('%02d', month(T)), '.', sprintf('%02d', T.Day)], ...
%                       [sprintf('%02d', T.Hour), sprintf('%02d', (T.Minute))]);
% if ~exist(OUT_FOLDER, 'dir')
%     mkdir(OUT_FOLDER);
% end

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% libigl MIQ (greedy rounding) %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    print_header('Solve with libigl MIQ')
    res1 = nrosy_wrapper(fp, f0, v0, N);
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %% Solve with Gurobi %%
    %%%%%%%%%%%%%%%%%%%%%%%
    print_header('Solve with Gurobi')
    tic
    %P = [d0'*d0, zeros(m.nV); zeros(m.nV), eye(m.nV)];
    P = [d0'*d0, zeros(m.nV); zeros(m.nV), zeros(m.nV)];
    q = zeros(2*m.nV, 1);
    r = 0;
    A = [d0'*d0, -(pi/2)*eye(m.nV);
         zeros(1, m.nV), ones(1, m.nV);
         ones(1, m.nV), zeros(1, m.nV)];
    b = [-Ad; 4*xi; 0];
    cvx_solver gurobi
    cvx_begin
    cvx_solver_settings('TimeLimit', 10)
        variable x_cvx(2*m.nV)
        variable z(m.nV) binary
        minimize (0.5*x_cvx'*P*x_cvx)
        subject to
            A*x_cvx == b
            x_cvx(m.nV+1:2*m.nV) == z
            x_cvx(m.nV+1:2*m.nV) <= 10
            x_cvx(m.nV+1:2*m.nV) >= -10
            %sum(x_cvx(m.nV+1:2*m.nV)) == N*xi
            %sum(x_cvx(1:m.nV)) == 0
            %ones(m.nV, 1)'*x_cvx(m.nV+1:2*m.nV) <= 10
            %ones(m.nV, 1)'*x_cvx(m.nV+1:2*m.nV) >= -10
    cvx_end
    toc
    
    res2.u = x_cvx(1:m.nV);
    res2.x = d0*res2.u;
    res2.k = x_cvx(m.nV+1:end);
    [res2.nrosy, res2.theta] = create_direction_field(m, res2.x, f0, theta0, N);
    res2.E = norm(res2.x)^2;
    res2.S_inds = find(res2.k > EPS);
    
    check_norm('d0''*res2.x-(pi/2)*res2.k','-Ad');
    
    res2.E
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run TCoDS with k from previous step %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    inds = res2.S_inds;
    I = [inds, res2.k(inds)/N];
    %theta0 = [res1.local_frames(1,:); res1.local_frames(m.nF+1,:)] * v0';
    %theta0 = atan2(theta0(2), theta0(1));
    %theta0 = 0;
    [res3.x, res3.E] = create_trivial_connection(m, I, N, VERBOSE);
    [res3.nrosy, res3.theta] = create_direction_field(m, res3.x, f0, theta0, N, res1.local_frames, res1.frame_diffs); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %% CoV from cvx to MIQ %%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    [res4.theta, res4.p] = TCoDS_to_MIQ(res2.x, res2.k, res1.frame_diffs, theta0, Ad, d0, d1, res1.period_is_fixed);
    res4.E = E_MIQ(m, res4.theta, res1.frame_diffs, res4.p, N);
    res4.nrosy = angles_to_nrosy(res4.theta, res1.local_frames, N);
    
    %%%%%%%%%%%%%%
    %% MIQ ADMM %%
    %%%%%%%%%%%%%%
    %      - k - | -- u --  
    % P = [  0   | d0'*d0 ]
    %     [  0   | 0      ]
    %
    %     [ -pi/2*I | d0'*d0 ]
    % A = [ 0       | 0      ]
    %     [ 1       | 0      ]
    %
    %     [-Ad   ]
    % b = [ 4*xi ]
    %     [ 0    ]
    xtmp = x_cvx;
    xtmp(1:m.nV) = x_cvx(m.nV+1:end);
    xtmp(m.nV+1:end) = x_cvx(1:m.nV);
    P_admm = ([zeros(m.nV), d0'*d0; 
               zeros(m.nV), zeros(m.nV)]);
    q_admm = zeros(2*m.nV, 1);
    A_admm = [-(pi/2)*eye(m.nV), d0'*d0;
              ones(1, m.nV), zeros(1, m.nV);
              zeros(1, m.nV), ones(1, m.nV)];
    r_admm = 0;
    b_admm = [-Ad; 4*xi; 0];
    Problem.objective.P = full(P_admm);
    Problem.objective.q = q_admm;
    Problem.objective.r = r_admm;
    %Problem.constraint.A = full(A);
    %Problem.constraint.b = b;
    Problem.constraint.A = full(A_admm);
    Problem.constraint.b = b_admm;
    Problem.constraint.l1 = m.nV;
    Problem.constraint.l2 = 0;
    Problem.constraint.l3 = 0;
    Problem.constraint.l4 = 0;
    
    maxiter = 5000; % max number of iterations
    repeat = 10; % number of initializations
    xs = solver_miqp_admm(Problem, 1, maxiter, repeat);

    f_admm = Inf;
    x_admm = zeros(2*m.nV,1);
    res_thr = 1e-3; % feasibilty threshold 
    thr = zeros(repeat*maxiter, 1);
    f = zeros(repeat*maxiter, 1);
    i = 1;
    normA = norm(full(A_admm));
    for j=1:repeat
        for k=1:maxiter
            x = xs(:,(j-1)*maxiter+k);
            thr(i) = norm(A_admm*x-b_admm);
            f(i) = 0.5*x'*P_admm*x + q_admm'*x + r_admm;
            %fprintf('%f, %f\n', thr(i), f(i))
            if thr(i) < res_thr
                %f = 0.5*x'*P*x + q'*x + r;
                if f(i) < f_admm
                    f_admm = f(i);
                    x_admm = x;
                end
            end
            i = i + 1;
        end
    end
    
    %assert(f_admm < Inf);
    min(thr)
    fprintf('The best value found by ADMM is %3s\n', f_admm);
    
    res5.u = x_admm(m.nV+1:end);
    res5.k = x_admm(1:m.nV);
    res5.x = d0*res5.u;
    
    check_norm('d0''*res5.x-(pi/2)*res5.k','-Ad')

    %assert(norm(res4.x - res2.x) < EPS)
    %norm(d1*res4.x)
    
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
    hold on;
    draw_constraints(m, f0, v0);
    draw_direction_field(m, res4.nrosy, res2.S_inds, scale)
    hold off
    title({'TCoDScvx(x,k) --> MIQ(p,theta)', ['E = ', num2str(res4.E)]})

    hlink = linkprop(ha, {'CameraPosition', 'CameraUpVector'});
    figi = figi + 1;
end