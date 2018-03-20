function [] = exp_2017_08_01_1234_globally_optimal()
ME = [];
try
%MESHES = {'sphere_s0.off'};
ABOUT = 'STAR6 Figure\r\n\r\n';

%MESHES = {'sphere_s0.off', ...
%    'bumpy.off', ...
%    'round_cuber.off', ...
%    'cow.off', ...
%    'bunny.off', ...
%    'phands.off', ...
%    'torus_fat_r2.off'};

%MESHES = {'sphere_s0.off'};
MESHES = {'squares5_row.off'};

%view_angles = {[-37.5,30], ...
%    [-65,50], ...
%    [40.9,-38.8], ...
%    [-179.1,-89.2], ...
%    [-161.5,10], ...
%    [95.7,20.4], ...
%    [18.5,-39.6]};

VERBOSE = true;
EPS = 1e-9;
LABEL_F = false;
LABEL_E = true;
LABEL_V = true;
DRAW_LOCAL_FRAMES = false;

theta0 = 0;

RUN_LIBIGL_MIQ = true;
PLOT = true;
SAVE = true;

if PLOT && SAVE
    OUT_FOLDER = create_time_stamped_folder(fullfile('..', 'results'));
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
    degree = 4;
    f0 = [1, 19]; % Starting face
    
    %v0 = mesh.V(mesh.F(1, 2), :) - mesh.V(mesh.F(1, 1), :);
    %MIQ = nrosy_wrapper(fp, [1], v0, degree);   
    
    %theta0 = [pi/4, pi/4+pi/2];
    theta0 = [pi/4, 0];
    %theta0 = [0.1, 0];
    v0 = [cos(theta0(1)), sin(theta0(1)), 0;
          cos(theta0(2)), sin(theta0(2)), 0];
    %local_frames = MIQ.local_frames;
    %v0 = angles_to_ffield(theta0, local_frames, 1);
    %theta0 = [pi/4, -pi/4-pi/2];
    %theta0 = [pi/4, pi/4+pi/2+0.1];
    %theta0 = [0, 0];

    log_and_print(LOG, 'degree : %d\r\n', degree);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % libigl MIQ (greedy rounding) %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if RUN_LIBIGL_MIQ
        log_and_print(LOG, '\r\nSolving with libigl MIQ...\r\n');
        MIQ = nrosy_wrapper(fp, f0, v0, degree);
        log_and_print(LOG, '%s\r\n', MIQ.result);
        log_and_print(LOG, 'Status: %d\r\n', MIQ.status);
    end
    
    %%%%%%
    % TC %
    %%%%%%
    log_and_print(LOG, '\r\nSolving with TCoDS...\r\n');
    [local_frames, ~] = create_local_frames(mesh);
    theta0 = gangles_to_local_angles(theta0, local_frames, 1, f0);
    %v0_local = angles_to_ffield(theta0, local_frames, 1, [1,19]);
    %theta0 = [atan2(v0_local(1, 2), v0_local(1, 1)), ...
    %          atan2(v0_local(2, 2), v0_local(2, 1))];
    J = [5, 7, 9, 14, 16, 21, 23, 28, 30, 35, 36, 33, 31, 26, 24, 19, 17, 12, 10, 3];
    I = ones(size(J));
    V = -ones(size(J));
    %V = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,1,1];
    boundary_constraints.A = sparse(I, J, V, 1, mesh.nE);
    boundary_constraints.b = [0];
    %boundary_constraints = [];
    TC = TCODS(...
        mesh, ...
        MIQ.S, ...
        f0, ...
        theta0, ...
        degree, ...
        'BoundaryConstraints', boundary_constraints, ...
        'Verbose', VERBOSE);
    
    %%%%%%%%%%%%%%%%%%%%
    % Globally optimal %
    %%%%%%%%%%%%%%%%%%%%
%     log_and_print(LOG, '\r\nSolving with GO...\r\n');
%     A = sparse(mesh.nE, mesh.nF);
%     EF = mesh.EFAdj;
%     [local_frames, diffs] = create_local_frames(mesh);
%     r = exp(1i*degree*diffs);
%     for eid = 1:mesh.nE
%         if mesh.isBoundaryEdge(eid)
%             %error('Boundary edges not implemented')
%             continue
%         end
%         
%         fi = EF(eid, 1);
%         fj = EF(eid, 2);
%         A(eid, fi) = r(eid);
%         A(eid, fj) = -1;
%     end
%       
%     %[V, ~] = eigs(A'*A, 1, 'sm');
%     [V, D] = eig(full(A'*A));
%     V = V(:, 1);
%     
%     GO.theta = angle(V) / degree;
%     GO.ffield = angles_to_ffield(GO.theta, local_frames, degree);
%     GO.degree = degree;
%     GO.p = zeros(mesh.nE, 1);
%     for eid = 1:mesh.nE
%         if mesh.isBoundaryEdge(eid)
%             %error('Boundary edges not implemented')
%             continue
%         end
%         fi = EF(eid, 1);
%         fj = EF(eid, 2);
%         GO.p(eid) = round(-2/pi * (GO.theta(fi) + diffs(eid) - GO.theta(fj)));
%     end
%     GO.E = E_MIQ(mesh, GO.theta, diffs, GO.p, degree);
%     [GO.x, GO.k] = MIQ_to_TCoDS(GO.theta, GO.p, diffs, MIQ.Kg, d0, d1);
%     inds = find(abs(GO.k) > EPS);
%     GO.S = [inds, GO.k(inds)];
   
    %%%%%%%%%
    % Plots %
    %%%%%%%%%
    
    if SAVE && PLOT
        log_and_print(LOG, '\r\nPlotting...\r\n');
        %MeshVis.wfigs('abc', mesh);
        %MIQ.degree = 1;
        %TC.degree = 1;
        MeshVis.wfigs(fname{:}, ...
            mesh, ...
            'OutFolder', OUT_FOLDER, ...
            'Titles', {['MIQ ', num2str(MIQ.E)], ...
                       ['TCoDS ', num2str(TC.E)]}, ...
            'Nrosy', {MIQ, TC}, ...
            'View', [0, 90], ...
            'FaceAlpha', 0, ...
            'Scale', 1/8, ...
            'LabelEdges', {true, TC.x}, ...
            'LabelVertices', LABEL_V, ...
            'LabelFaces', LABEL_F, ...
            'LocalFrames', DRAW_LOCAL_FRAMES);
    elseif PLOT
        figure
        subplot(221)
        MeshVis.plot(mesh, 'nrosy', MIQ);
        title(['MIQ ', num2str(MIQ.E)])
        subplot(222)
        MeshVis.plot(mesh, 'nrosy', TC);
        title(['TCoDS ', num2str(TC.E)])
        subplot(223)
        %MeshVis.plot(mesh, 'nrosy', GO);
        %title(['GO ', num2str(GO.E)])
    end
    
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

end