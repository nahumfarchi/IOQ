% 3x3 grid
%   (1,pi/4) -> (2,pi/4)   -> (3,pi/4)
%       |           |           |
%   (4,pi/4) -> (5,-5pi/6) -> (6,pi/4)
%       |           |           |
%   (7,pi/4) -> (8,pi/4)   -> (9,pi/4)
%
% Solution should be (with N=1):
%   n_5 = 1
%   all other variables = 0
%
% See scripts\PI_grid_3x3_example.PNG


%%%%%%%%%
% Setup %
%%%%%%%%%
N = 1; % Field degree

F = [1, 4, 5, 2;
     2, 5, 6, 3;
     4, 7, 8, 5;
     5, 8, 9, 6];
V = [0, 0, 0;
     0, 1, 0;
     0, 2, 0;
     1, 0, 0;
     1, 1, 0;
     1, 2, 0;
     2, 0, 0;
     2, 1, 0;
     2, 2, 0];
E = [1, 2;
     1, 4;
     2, 3;
     2, 5;
     3, 6;
     4, 5;
     4, 7;
     5, 6;
     5, 8;
     6, 9;
     7, 8;
     8, 9];
nF = size(F, 1);
nV = size(V, 1);
nE = size(E, 1);

wrapped = (pi/4)*ones(nV, 1);
wrapped(5) = -5*pi/6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Edgelist phase unwrapping %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[unwrapped, n, P, Q] = unwrap_phases_grid(V, E, wrapped, N);
assert(n(5)-1 < 1e-8)
assert(~any(n(1:4) > 1e-8))
assert(~any(n(6:end) > 1e-8))


%%%%%%%%%%%%%%%%%%
% Costantini DEC %
%%%%%%%%%%%%%%%%%%

% Edge orientation is (i,j), i<j
% Face orientation is CCW
d_0 = zeros(nE, nV);
for e = 1:nE
    v1 = E(e, 1);
    v2 = E(e, 2);
    d_0(e, v2) = 1;
    d_0(e, v1) = -1;
end
% Edge  1   2   3   4   5   6   7   8   9   10   11   12   |
%       ---------------------------------------------------| F
d_1 = [-1,  1,  0, -1,  0,  1,  0,  0,  0,  0,   0,   0;   % 1
        0,  0, -1,  1, -1,  0,  0,  1,  0,  0,   0,   0;   % 2
        0,  0,  0,  0,  0, -1,  1,  0, -1,  0,   1,   0;   % 3
        0,  0,  0,  0,  0,  0,  0, -1,  1, -1,   0,   1];  % 4
    
wrapped_tilde = d_0*wrapped;

% find integers such that wrapped_tilde \in [-pi,pi)
tmp = (wrapped_tilde+sign(wrapped_tilde)*pi) / (2*pi);
n_tilde = zeros(size(tmp));
n_tilde(tmp>0) = -floor(tmp(tmp>0));
n_tilde(tmp<0) = -ceil (tmp(tmp<0));

wrapped_tilde = wrapped_tilde + 2*pi*n_tilde;
assert(~(any(wrapped_tilde < -pi | wrapped_tilde > pi)))

% Solve min |R|_1
% s.t.  d_1(R) = -d_1*n_tilde
%       R is integer
cost = ones(2*nE, 1);
lb = zeros(2*nE, 1);
ub = 50*ones(2*nE, 1);
A = [d_1, -d_1]; % R = x^+ - x^-
b = -d_1*n_tilde;
[x, fval, exitflag] = linprog(cost, [], [], A, b, lb, ub);

% Since d_1*n_tilde=0 we get x=0. So the gradient is just wrapped_tilde.
% Looking at this gradient, we can see that the only vertex that changes 
% is v5:
%   v5 = v1+e1+e4 = v1+e4
unwrapped_costantini = wrapped;
unwrapped_costantini(5) = unwrapped_costantini(1)+wrapped_tilde(4);
% and then we get the same solution
assert(norm(unwrapped-unwrapped_costantini) < 1e-8)












%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cost = ones(nE, 1);
% ub = 50*ones(nE, 1);
% lb = -50*ones(nE, 1);
% A = d_1;
% b = -d_1*n_tilde;
% %[x, fval, exitflag] = linprog(cost, [], [], A, b, lb, ub);
% [R, fval, exitflag] = intlinprog(cost, 1:nE, [], [], A, b, lb, ub);


% d_11 = zeros(nF, nE);
% for f = 1:nF
%     f_verts = F(f, :);
%     n_verts = length(f_verts);
%     for i = 1:n_verts
%         v1 = f_verts(i);
%         v2 = f_verts(1+mod(i, n_verts));
%         [C, IA, IB] = intersect(E, [v1, v2], 'rows');
%         if ~isempty(C)
%             d_11(f, IA) = 1;
%         else
%             [C, IA, IB] = intersect(E, [v2, v1], 'rows');
%             d_11(f, IA) = -1;
%         end
%     end
% end
