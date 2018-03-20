% Solve argmin_{u,k} \|d_0u\|_2^2
%       s.t. d_0^Td_0u = -A_d + (pi/2)k
%            k \in \mathbb{Z}
%            \sum k(i) = \sum A_d(i)
%            \sum u(i) = 0

%% Settings
%FNAME = 'sphere_s0.off';
%FNAME = 'sphere_s1.off';
%FNAME = 'ellipsoid_s4r.off';
%FNAME = 'round_cuber.off';
%FNAME = 'sphere_s0.off';
%FNAME = 'torus_fat_r2.off';
%FNAME = 'bunny.off';

%MESHES = {'sphere_s0.off', 'ellipsoid_s4r.off', 'round_cuber.off', 'torus_fat_r2.off', 'bunny.off', 'phands.off'};
MESHES = {'sphere_s0.off'};

VERBOSE = true;
EPS = 1e-9;
EDGE_ALPHA = 0.2;
RESOLUTION = 1024;
AR = 1;

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
    
    %% Minimize using CoMISo
    % Constraints
    C = [d0'*d0, -(pi/2)*eye(m.nV); ...
        zeros(1,m.nV), ones(1,m.nV); ...
        ones(1,m.nV), zeros(1,m.nV)];
    bc = [Ad; -N*xi; 0];
    C = [C, bc]; % This is how CoMISo wants it (I think)
    
    % Remove half of the constraints
    %C = [C(1:floor(m.nV/2),:); ...
    %    C(end-1,:); ...
    %    C(end,:)];
    
    % Leave only 10 constraints
    %C = [C(1:10,:); ...
    %    C(end-1:end,:)];
    
    % Gradient matrix
    %A = [2*(d0'*d0)*(d0'*d0), zeros(m.nV, m.nV); ...
    %     zeros(m.nV, m.nV), 2*eye(m.nV)];
    
    %A = [d0'*d0, zeros(m.nV, m.nV);
    %     zeros(m.nV, m.nV), eye(m.nV)];
         %zeros(m.nV, m.nV + m.nV)];
    A = [d0'*d0, zeros(m.nV);
         zeros(m.nV), zeros(m.nV)];
    B = [2*(d0'*d0)*(d0'*d0), zeros(m.nV);
         zeros(m.nV), 2*zeros(m.nV)];
    ba = zeros(m.nV + m.nV, 1);
    bb = zeros(m.nV + m.nV, 1);
    % Integer variables
    idx_to_round = m.nV+1 : 2*m.nV;
    %idx_to_round = [2*m.nV];
    
    %% Check grad
    f = @(x) E(x, A, ba, B, bb);
    fcon = @(x) deal([], C(:,1:end-1)*x-C(:,end));
    %check_grad(@(x) ftmp(x), fcon, 4);
    %check_grad(f, fcon, size(A,2));
    %check_grad(@(x) norm(A*x-ba)^2, @(x) (B*x - bb), @(x) [], @(x) C(:,1:end-1)*x-C(:,end), size(A,2));
    %check_grad(@(x) E_2017_07_16_1105(x, d0, Ad), fcon, size(A,2));
    check_grad(@(x) ftmp(x, d0, Ad), fcon, size(A,2));
    
    %% Solve
    x = CoMISo_wrapper(full(C), full(B), ba, idx_to_round, true); % TODO add option for sparse matrix input to wrapper function
    u = x(1:m.nV);
    k = x(m.nV+1:end);
    connection = d0*u;
    
    norm(d0'*connection - (-Ad+(pi/2)*k))
    sum(k)
    
    theta0 = 0;
    [nrosy, theta] = create_direction_field(m, connection, f0, theta0, N);
    
    %% Run TCoDS with k from previous step
    inds = find(k > EPS);
    I = [inds, k(inds)/N];
    %theta0 = [res1.local_frames(1,:); res1.local_frames(m.nF+1,:)] * v0';
    %theta0 = atan2(theta0(2), theta0(1));
    %theta0 = 0;
    [res2.x, res2.E] = create_trivial_connection(m, I, N, VERBOSE);
    [res2.nrosy, res2.theta] = create_direction_field(m, res2.x, f0, theta0, N, res1.local_frames, res1.frame_diffs); %TODO use same angles as in MIQ above

    %assert(norm(res4.x - res2.x) < EPS)
    %norm(d1*res4.x)
    
    %% Plot
    figure
    draw_direction_field(m, nrosy, [], 0.5*m.avg_length);
%     
%     %% libigl MIQ (greedy rounding)
%     print_header('libigl MIQ')
%     res1 = nrosy_wrapper(fp, f0, v0, N);
%     
%     %%
%     func = (d1*wrapToPi(res1.frame_diffs) + (pi/2)*d1*mod(res1.p,4));
% 
%     figure
%     m.draw(func, 'EdgeAlpha', EDGE_ALPHA);
%     title('d_1r+(\pi/2)d_1p')
%     colorbar
%     
%     [~, name_only, ~] = fileparts(fname{:});
%     pngname = fullfile(OUT_FOLDER, [name_only, '.png']);
%     %pngname = sprintf('%s_%04d.png',filename,i);
%                 
%     dpi = get(0, 'ScreenPixelsPerInch');
%     in = RESOLUTION/dpi;
% 
%     fig = gcf;
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 AR*in in];
%     fig.PaperPositionMode = 'manual';
%     print(pngname, '-dpng', '-r0')
end
