%% =========================================================================
%  lamg_example.m
%
%  LAMG Example usage: solve a graph Laplacian system with a random
%  compatible right-hand side.
%  =========================================================================

fprintf('Setting up problem\n');
%SIZE = 10000;
%g = Graphs.grid('fd', [SIZE SIZE], 'normalized', true);
%A = g.laplacian;                % Zero row-sum (Neumann B.C.)
%b = rand(size(A,1), 1);
%b = b - mean(b);                % Make RHS compatible with A's null space
mesh = Mesh('../data/ashish_nob/armadillo.off');

d0 = get_exterior_derivatives(mesh);

p = colamd(d0);
pt(p) = 1:length(p);
tColamd = timeit(@() colamd(d0));
d0p = d0(:,p);
A = d0'*d0;
Ap = d0p'*d0p;


% A = d0'*d0;
% p = colamd(A);
% tColamd = timeit(@() colamd(A));
% pt(p) = 1:length(p);
% Ap = A(:, p);

b = rand(size(A,1), 1);
b = b - mean(b);

%A = d0'*d0;
%p = symamd(A);
%tColamd = timeit(@() symamd(A));
%pt(p) = 1:length(p);
%Ap = A(p, p);


inputType = 'laplacian';        % The input matrix A is a graph Laplacian
solver = 'lamg';                % Or 'cmg', or 'direct'

setup_time = [];
solve_time = [];
total_time = [];
err = [];
names = {};

%% ---------------------------------------------------------------------
%  Setup phase: construct a LAMG multi-level hierarchy (A)
%  ---------------------------------------------------------------------
fprintf('Setting up solver %s\n', solver);
lamg    = Solvers.newSolver(solver, 'randomSeed', 1);
%tStart  = tic;
setup   = lamg.setup(inputType, A);
%tSetup  = toc(tStart);
tSetup = timeit(@() lamg.setup(inputType, A));
setup_time(end+1) = tSetup;

%% ---------------------------------------------------------------------
%  Solve phase: set up a random compatible RHS b and solve A*x=b. You can
%  repeat this phase for multiple b's without rerunning the setup phase.
%  ---------------------------------------------------------------------
setRandomSeed(now);
% Turn on debugging printouts during the run
core.logging.Logger.setLevel('lin.api.AcfComputer', core.logging.LogLevel.DEBUG);
fprintf('Solving A*x=b\n');
%tStart = tic;
[x, ~, ~, details] = lamg.solve(setup, b, 'errorReductionTol', 1e-4);
%tSolve = toc(tStart);
tSolve = timeit(@() lamg.solve(setup, b, 'errorReductionTol', 1e-4));

solve_time(end+1) = tSolve;
total_time(end+1) = tSetup + tSolve;
err(end+1) = norm(A*x - b) / norm(b);

names{end+1} = 'LAMG (A)';
% Turn printouts off
core.logging.Logger.setLevel('lin.api.AcfComputer', core.logging.LogLevel.INFO);
fprintf('------------------------------------------------------------------------\n');

%% ---------------------------------------------------------------------
%  Solve with ichol (A)
%  ---------------------------------------------------------------------
disp('ichol')
%tStart = tic;
Aichol = ichol(A,struct('type','ict','droptol',1e-4,'michol','off'));
%tSetupIchol = toc(tStart);
tSetupIchol = timeit(@() ichol(A,struct('type','ict','droptol',1e-4,'michol','off')));

%tStart = tic;
%xichol = Achol' \ (Achol \ b);
%tSolveIchol = toc(tStart);
tSolveIchol = timeit(@() Aichol' \ (Aichol \ b));
xichol = Aichol' \ (Aichol \ b);


setup_time(end+1) = tSetupIchol;
solve_time(end+1) = tSolveIchol;
total_time(end+1) = tSetupIchol + tSolveIchol;
err(end+1) = norm(A*xichol - b) / norm(b);
names{end+1} = 'ichol (A)';

%% ---------------------------------------------------------------------
%  Solve with chol (A)
%  ---------------------------------------------------------------------
disp('chol')
tStart = tic;
[Achol, ~] = chol(A);
tSetupChol = toc(tStart);
%tSetupChol = timeit(@() chol(A));

%tStart = tic;
%xichol = Achol' \ (Achol \ b);
%tSolveIchol = toc(tStart);
tSolveChol = timeit(@() Achol \ (Achol' \ b));
xchol = Achol \ (Achol' \ b);

setup_time(end+1) = tSetupChol;
solve_time(end+1) = tSolveChol;
total_time(end+1) = tSetupChol + tSolveChol;
err(end+1) = norm(A*xchol - b) / norm(b);
names{end+1} = 'chol (A)';


%% 
%---------------------------------------------------------------------
% Solve with CMG (A)
%---------------------------------------------------------------------
disp('cmg')
tol = 1e-4;
numIterations = 300;

%tStart = tic;
pfun = cmg_sdd(A);
%tSetupCMG = toc(tStart);
tSetupCMG = timeit(@() cmg_sdd(A));

%tStart = tic;
%x = pcg(A, b, tol, numIterations, pfun);
%tSolveCMG = toc(tStart);
tSolveCMG = timeit(@() pcg(A, b, tol, numIterations, pfun));
xpcg = pcg(A, b, tol ,numIterations, pfun);

setup_time(end+1) = tSetupCMG;
solve_time(end+1) = tSolveCMG;
total_time(end+1) = tSetupCMG + tSolveCMG;
err(end+1) = norm(A*xpcg - b) / norm(b);
names{end+1} = 'CMG+PCG (A)';

%%
%---------------------------------------------------------------------
% Setup phase: construct a LAMG multi-level hierarchy (Ap)
%---------------------------------------------------------------------
fprintf('Setting up solver %s\n', solver);
lamg    = Solvers.newSolver(solver, 'randomSeed', 1);
%tStart  = tic;
setup   = lamg.setup(inputType, Ap);
%tSetup  = toc(tStart);
tSetup = timeit(@() lamg.setup(inputType, Ap)) + tColamd;
setup_time(end+1) = tSetup;

%%
%---------------------------------------------------------------------
% Solve phase: set up a random compatible RHS b and solve A*x=b. You can
% repeat this phase for multiple b's without rerunning the setup phase.
%---------------------------------------------------------------------
setRandomSeed(now);
% Turn on debugging printouts during the run
core.logging.Logger.setLevel('lin.api.AcfComputer', core.logging.LogLevel.DEBUG);
fprintf('Solving A*x=b\n');
%tStart = tic;
[x, ~, ~, details] = lamg.solve(setup, b(p), 'errorReductionTol', 1e-4);
x = x(pt);
%tSolve = toc(tStart);
tSolve = timeit(@() lamg.solve(setup, b(p), 'errorReductionTol', 1e-4));

solve_time(end+1) = tSolve;
total_time(end+1) = tSetup + tSolve;
err(end+1) = norm(A*x - b) / norm(b);

names{end+1} = 'LAMG (Ap)';
% Turn printouts off
core.logging.Logger.setLevel('lin.api.AcfComputer', core.logging.LogLevel.INFO);
fprintf('------------------------------------------------------------------------\n');

%% ---------------------------------------------------------------------
%  Solve with ichol (Ap)
%  ---------------------------------------------------------------------
disp('ichol')
%tStart = tic;
Aichol = ichol(Ap,struct('type','ict','droptol',1e-4,'michol','off'));
%tSetupIchol = toc(tStart);
tSetupIchol = timeit(@() ichol(Ap,struct('type','ict','droptol',1e-4,'michol','off'))) + tColamd;

%tStart = tic;
%xichol = Achol' \ (Achol \ b);
%tSolveIchol = toc(tStart);
currentParms = spparms;
spparms('autoamd',0);
spparms('autommd',0);

tSolveIchol = timeit(@() Aichol' \ (Aichol \ b(p)));
z = Aichol' \ (Aichol \ b(p));
xichol = z(pt);

spparms(currentParms);

setup_time(end+1) = tSetupIchol;
solve_time(end+1) = tSolveIchol;
total_time(end+1) = tSetupIchol + tSolveIchol;
%err(end+1) = norm(Ap*xichol - b) / norm(b);
err(end+1) = norm(A*xichol - b) / norm(b);
names{end+1} = 'ichol (Ap)';

%% ---------------------------------------------------------------------
%  Solve with chol (Ap)
%  ---------------------------------------------------------------------
disp('chol')
tStart = tic;
%[Achol, ~] = chol(Ap);
[Achol, ~, S] = chol(A);
tSetupChol = toc(tStart);
%tSetupChol = timeit(@() chol(Ap)) + tColamd;

%tStart = tic;
%xichol = Achol' \ (Achol \ b);
%tSolveIchol = toc(tStart);
c = S'*b;
tSolveChol = timeit(@() Achol \ (Achol' \ c));
z = Achol \ (Achol' \ c);
xchol = S * z;

setup_time(end+1) = tSetupChol;
solve_time(end+1) = tSolveChol;
total_time(end+1) = tSetupChol + tSolveChol;
err(end+1) = norm(A*xchol - b) / norm(b);
names{end+1} = 'chol (Ap)';


%% ---------------------------------------------------------------------
%  Solve with CMG (Ap)
%  ---------------------------------------------------------------------
disp('cmg')
tol = 1e-4;
numIterations = 300;

%tStart = tic;
pfun = cmg_sdd(Ap);
%tSetupCMG = toc(tStart);
tSetupCMG = timeit(@() cmg_sdd(Ap)) + tColamd;

%tStart = tic;
%x = pcg(A, b, tol, numIterations, pfun);
%tSolveCMG = toc(tStart);
tSolveCMG = timeit(@() pcg(Ap, b(p), tol, numIterations, pfun));
xpcg = pcg(Ap, b(p), tol ,numIterations, pfun);
xpcg = xpcg(pt);

setup_time(end+1) = tSetupCMG;
solve_time(end+1) = tSolveCMG;
total_time(end+1) = tSetupCMG + tSolveCMG;
err(end+1) = norm(A*xpcg - b) / norm(b);
names{end+1} = 'CMG+PCG (Ap)';

%% Table
T = table(setup_time', solve_time', total_time', err', 'RowNames', names, 'VariableNames', {'SetupTime', 'SolveTime', 'TotalTime', 'Error'})

return

%% 
%---------------------------------------------------------------------
% Solve with lsqminnorm
%---------------------------------------------------------------------
%tLsqminnorm = timeit(@() lsqminnorm(A, b));

%% 
%---------------------------------------------------------------------
% Solve with colamd
%---------------------------------------------------------------------
A = d0'*d0;
tic
p = colamd(A);
toc
figure;
subplot(121); spy(A); title('A')
subplot(122); spy(A(:,p)); title('A(:,p)')

[Achol1, ~] = chol(A);
[Achol2, ~] = chol(d0(:,p)'*d0(:,p));
figure
subplot(121); spy(Achol1); title('Achol1')
subplot(122); spy(Achol2); title('Achol2')

Aichol1 = ichol(A, struct('type','ict','droptol',1e-4,'michol','off'));
Aichol2 = ichol(d0(:,p)'*d0(:,p), struct('type','ict','droptol',1e-4,'michol','off'));
figure
subplot(121); spy(Aichol1); title('Aichol1')
subplot(122); spy(Aichol2); title('Aichol2')

%%
%---------------------------------------------------------------------
% Display statistics
%---------------------------------------------------------------------
disp(setup);
tMvm    = mvmTime(A, 5);
nnz     = numel(nonzeros(A));

fprintf('\n');
fprintf('MVM time [sec]       elapsed %.3f, per nonzero %.2e\n', ...
    tMvm, tMvm / nnz);
fprintf('Setup time [sec]     elapsed %.3f, per nonzero %.2e, in MVM %.2f\n', ...
    tSetup, tSetup / nnz, tSetup / tMvm);
fprintf('Solve time [sec]     elapsed %.3f, per nonzero %.2e, in MVM %.2f\n', ...
    tSolve, normalizedSolveTime(tSolve, details.errorNormHistory) / nnz, tSolve / tMvm);
fprintf('total                elapsed %.3f\n', tSetup+tSolve);
fprintf('|A*x-b|/|b|          %.2e\n', norm(A*x-b)/norm(b));
if (isfield(details, 'acf'))
    fprintf('Convergence factor   %.3f\n', details.acf);
end
fprintf('\n');
fprintf('ichol setup time     elapsed %.3f\n', tSetupIchol);
fprintf('ichol solve time     elapsed %.3f\n', tSolveIchol);
fprintf('total                elapsed %.3f\n', tSetupIchol+tSolveIchol);
fprintf('|A*x-b|/|b|          %.2e\n', norm(A*xichol-b)/norm(b));
fprintf('\n');
fprintf('CMG setup time       elapsed %.3f\n', tSetupCMG);
fprintf('CMG solve time       elapsed %.3f\n', tSolveCMG);
fprintf('total                elapsed %.3f\n', tSetupCMG+tSolveCMG);
fprintf('|A*x-b|/|b|          %.2e\n', norm(A*xpcg-b)/norm(b));
fprintf('\n');
