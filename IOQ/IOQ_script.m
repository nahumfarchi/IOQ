%%%%%
FACE0 = 1;        % starting face for tcods
THETA0 = nan;     % (use global constraint instead)
GVEC = [1, 0, 0]; % constraint vector (in global coordinates)
DEGREE = 4;       % works only for degree 4 atm

disp('Loading mesh...')
mesh = Mesh('../../data/bunny.off');
V = mesh.V;
F = mesh.F;
% or: 
%mesh = Mesh(V, F);

nv = mesh.nV; nf = mesh.nF;


%%%%%
disp('Running IOQ...')
k = IOQ_cpu(V, F, 'Laplacian', 'cot', 'Plot', true);


%%%%%
disp('Creating direction field...')
[A, K, d0, d1, H] = tcods_gsystem(V, F);
n = size(A, 2);
b = K - (pi / 2) * k;
x = lsqlin(speye(n, n), zeros(n, 1), [], [], A, -b, -inf(n, 1), inf(n, 1));

assert(norm(A*x - (-b)) < 1e-10)
assert(norm(d1*x) < 1e-10)
%assert(abs(rank(full([A; d1])) - min(size([A; d1]))) < 1e-10)

fprintf('E : %g\n', norm(x)^2);
[ffield] = connection_to_nrosy(...
                mesh, ...
                x, ...
                FACE0, ...
                THETA0, ...
                DEGREE, ...
                'gConstraintVec', GVEC);


%%%%%
disp('Plotting...')
mesh.ffield_vectors = ffield;
mesh.degree = DEGREE;
mesh.vert_sing = [find(k), k(k~=0)];
figure()
mesh.draw()
