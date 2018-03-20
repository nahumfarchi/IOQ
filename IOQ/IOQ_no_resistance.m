function [alpha, beta, elapsed_total, E] = IOQ_no_resistance(verts, faces, varargin)

parser = inputParser;

addOptional(parser, 'Iterations', 1000);
addOptional(parser, 'NSingularities', []);
addOptional(parser, 'Plot', false);
addOptional(parser, 'Tol', 1e-6);
addOptional(parser, 'beta_P', []);
addOptional(parser, 'Verbose', true);
addOptional(parser, 'Mesh', []);

parse(parser, varargin{:});
opt = parser.Results;
if isempty(opt.Mesh)
    mesh = Mesh(verts, faces);
else
    mesh = opt.Mesh;
end
nv = mesh.nV; ne = mesh.nE; nf = mesh.nF;

tic

d0 = get_exterior_derivatives(mesh);
L = d0'*d0;
F = factorize(L);
S = inverse(L);

alpha_g = (2/pi) * get_gaussian_curvature(mesh);
genus = round(1 - sum(alpha_g)/(4*pi));
ng2 = 2*genus;
if isempty(opt.NSingularities)
    c = round(abs(sum((alpha_g)))); % = 4 xi
else
    c = opt.NSingularities;
end

alpha = zeros(nv, 1);
n_pos_sing = (c + round(sum(alpha_g))) / 2;
n_neg_sing = (c - round(sum(alpha_g))) / 2;
inds_pos = randperm(nv, n_pos_sing);
inds_neg = randperm(nv, n_neg_sing);
alpha(inds_pos) = 1;
alpha(inds_neg) = -1;

if genus == 0
    beta = [];
else
    H = -mesh.H;
    z = wrapToPi(generator_angle_defects(mesh));
    K = [alpha_g; z];
    beta_g = K(nv+1:end);
    B = H - d0 * (F \ (d0' * H));
    A2 = (H' * d0) / F;
    y0 = (2/pi)*(A2*alpha_g - beta_g);
    beta = round(A2*alpha - y0);
end

% JLFac = 24;
% eps = 0.5;
% k = round(JLFac * log(nv) / eps^2);
% Y = 2*(rand(k, ne) > 0.5) - 1;
% Y = (1/sqrt(k)) * Y;
% Y = Y * d0;
% Ztilde = (F \ Y');
% EV = mesh.EVAdj;
% Rtilde = sparse(nv);
% for eid = 1:ne
%     v1 = EV(eid, 1);
%     v2 = EV(eid, 2);
%     rs = norm(Ztilde(v1, :) - Ztilde(v2, :))^2;
%     Rtilde(v1, v2) = rs;
%     Rtilde(v2, v1) = rs;
% end
%[Rtilde, ~] = resistance_distance(mesh, 0, 24, false);
Rtilde = sparse(nv, nv);

u = 2 * (F \ (alpha - alpha_g));
big_diag = sparse(1:nv, 1:nv, inf(nv, 1), nv, nv, nv);
E = [];
for iter = 1:opt.Iterations
    %E(iter) = (alpha-alpha_g)' * u;
    
    [min_vals, is] = min(bsxfun(@plus, u, -u') + big_diag);
    [~, inds] = sort(min_vals);
    js = inds(1:10);
    is = is(js);
    m_best = inf;
    i_best = nan;
    j_best = nan;
    for l = 1:10
        i = is(l);
        j = js(l);
        if Rtilde(i, j) == 0
            rs = S(i, i) + S(i, j) - 2*S(i, j);
            Rtilde(i, j) = rs;
            Rtilde(j, i) = rs;
        end
        m = Rtilde(i, j) + min_vals(j);
        if m < m_best
            m_best = m;
            i_best = i;
            j_best = j;
        end
    end
    
    %[min_val, i] = min(bsxfun(@plus, u, -u') + big_diag + Rtilde);
    %[min_val, j] = min(min_val);
    %i = i(j);
    alpha(i_best) = alpha(i_best) + 1;
    alpha(j_best) = alpha(j_best) - 1;
    if abs(m_best) < opt.Tol
        disp('Minimization complete because optimality tolerance was reached.')
        disp(['m : ', num2str(m_best)])
        break
    end
    
    m_best_hist(iter) = m_best;
    if iter > 2 && abs(m_best_hist(iter) - m_best_hist(iter-2)) < opt.Tol
        disp('Cycle detected.')
        break
    end
    
    u = 2 * (F \ (alpha - alpha_g));
end

elapsed_total = toc;

if iter == opt.Iterations
    fprintf('Optimization terminated because the iterations budget %d was reached\n', opt.Iterations);
    fprintf('min_val = %g\n', m_best);
end
    

end

