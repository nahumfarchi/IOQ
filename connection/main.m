%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trivial connections on discrete surfaces %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

degree = 1;                          % Field degree
S = [1, 1; 10, 1];                   % nSx2 Matrix of singularity pairs (face_id, ki)
CONSTRAINED_FACES = [5, 50, 100];    % A list of constrained face ids
CONSTRAINT_ANGLES = [pi, -pi, pi/2]; % A list of constraint angles

%MESHES = {'sphere_s0.off', 'ellipsoid_s4r.off', 'round_cuber.off', 'torus_fat_r2.off', 'bunny.off', 'phands.off'};
%MESHES = {'sphere_s0.off', 'bumpy.off', 'bunny.off', 'torus_fat_r2.off', 'phands.off'};
MESHES = {'sphere_s0.off'};
%MESHES = {'sphere_s0.off', ...
%    'bumpy.off', ...
%    'round_cuber.off', ...
%    'cow.off', ...
%    'bunny.off', ...
%    'phands.off', ...
%    'torus_fat_r2.off'};

view_angles = {[-37.5,30], ...
    [-65,50], ...
    [40.9,-38.8], ...
    [-179.1,-89.2], ...
    [-161.5,10], ...
    [95.7,20.4], ...
    [18.5,-39.6]};

VERBOSE = true;

if exist(fullfile('.', 'data'), 'dir')
    DATA_FOLDER = fullfile('.', 'data');
elseif exist(fullfile('..', 'data'), 'dir')
    DATA_FOLDER = fullfile('..', 'data');
else
    error('Could not find data folder.')
end
RESULTS_FOLDER = fullfile('.', 'results');
if ~exist(RESULTS_FOLDER, 'dir')
    mkdir(RESULTS_FOLDER);
end

k = 1;
T = zeros(length(MESHES), 2);
for fname = MESHES
    disp(fname)
    fp = fullfile(DATA_FOLDER, fname{:});
    
    mesh = Mesh();
    mesh.loadTM(fp);
    
    tic
    nrosy = TCODS(mesh, S, CONSTRAINED_FACES, CONSTRAINT_ANGLES, degree, 'Verbose', VERBOSE);
    T(k, 2) = toc;
    T(k, 1) = mesh.nV;
    
    out_fp = fullfile(RESULTS_FOLDER, fname{:});
    figure
    MeshVis.wfigs(out_fp, ...
        mesh, ...
        'Nrosy', {nrosy}, ...
        'OpenGL', true, ...
        'Scale', 1, ...
        'Montage', true, ...
        'View', view_angles{k}, ...
        'ConstrainedFaces', CONSTRAINED_FACES, ...
        'ConstraintAngles', CONSTRAINT_ANGLES);
    k = k + 1;
end

% Plot timings
figure
D = T(:, 1);
F = T(:, 2);
plot(D,F,'xr','MarkerSize', 15)
p = polyfit(D,F,1);
f = polyval(p,D);
hold on
plot(D,f,'k')
hold off
xlabel('|V|')
ylabel('Time (s)')