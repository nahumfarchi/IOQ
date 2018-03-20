%% ========================================================================
%  Visualize the resistance distance and the approximated resistance
%  distance on some mesh (relative to some vertex).
%  ========================================================================
plotting_defaults;
out_folder = 'results';
mkdir(out_folder);
EPS = 0.5;
JLFAC = 24;
USE_GPU = false;
SEED = 112;
SAVE = false;
VERT = 135;
%fp = '../../../data/bunny.off';
%fp = '../../../data/bunnies/bunny_26k_faces.off';
%fp = '../../../data/3holes_lo.off';
fp = '../../../data/elephant_r.off';
[~, meshname, ~] = fileparts(fp);
m = Mesh(fp); ne = m.nE; nv = m.nV;

if USE_GPU
    gd = gpuDevice();
else
    gd = [];
end

tic
Rtilde = resistance_distance(m, EPS, JLFAC, USE_GPU); 
if ~isempty(gd), wait(gd); end; elapsed_Rtilde = toc
tic
R = resistance_distance(m, 0);
if ~isempty(gd), wait(gd); end; elapsed_R = toc
fprintf('|R - Rtilde| = %.4g\n', norm(R - squareform(Rtilde), 'fro') / norm(R, 'fro'));

%
d0 = get_exterior_derivatives(m);
L = d0' * d0;
F = factorize(L);

alpha_g = (pi/2)*get_gaussian_curvature(m);
g = m.genus;
c = round(abs(sum(alpha_g)));

n_pos_sing = (c + round(sum(alpha_g))) / 2;
n_neg_sing = (c - round(sum(alpha_g))) / 2;
inds_pos = randperm(nv, n_pos_sing);
inds_neg = randperm(nv, n_neg_sing);
alpha = zeros(nv , 1);
alpha(inds_pos) = 1;
alpha(inds_neg) = -1;

%alpha = alpha_p;

u = 2 * (F \ (alpha - alpha_g));

%%
out_folder2 = fullfile(out_folder, 'schemaball');
mkdir(out_folder2);
%normalize_cnst = max(max(max(abs(-R+U))), max(max(abs(squareform(-Rtilde)+U))));

n_pts = 10;
pts = [VERT, setdiff(randperm(nv, n_pts), VERT)];
%pts = 1:nv;
dist = graphshortestpath(m.VVAdj, VERT);
dist = dist(pts);
[dist, inds] = sort(dist);

%%
U = bsxfun(@plus, u, -u');
X1 = R(pts, pts) + U(pts, pts);
X1 = X1 ./ max(abs(X1(:)));
X1 = X1(inds, inds);
h1 = schemaball(X1);
title({'Exact Resistance Distance', 'U+R'})
if SAVE
    filename = fullfile(out_folder2, 'U_R_00');
    export_fig([filename '.png']);
    export_fig([filename '.pdf']);
    saveas(gcf, filename, 'fig')
end

%U = bsxfun(@plus, u, -u');
X2 = squareform(Rtilde);
X2 = X2(pts, pts) + U(pts, pts);
X2 = X2 ./ max(abs(X2(:)));
X2 = X2(inds, inds);
h2 = schemaball(X2);
title({'Approximate Resistance Distance', '$U+\tilde R$'})
if SAVE
    filename = fullfile(out_folder2, 'U_R_01');
    export_fig([filename '.png']);
    export_fig([filename '.pdf']);
    saveas(gcf, filename, 'fig')
end

%%
X1 = R(pts, pts);
X1 = X1 ./ max(abs(X1(:)));
X1 = X1(inds, inds);
h1 = schemaball(X1);
title({'Exact Resistance Distance', 'R'})
if SAVE
    filename = fullfile(out_folder2, 'R_00');
    export_fig([filename '.png']);
    export_fig([filename '.pdf']);
    saveas(gcf, filename, 'fig')
end

%U = bsxfun(@plus, u, -u');
X2 = squareform(Rtilde);
X2 = X2(pts, pts);
X2 = X2 ./ max(abs(X2(:)));
X2 = X2(inds, inds);
h2 = schemaball(X2);
title({'Approximate Resistance Distance', '$\tilde R$'})
if SAVE
    filename = fullfile(out_folder2, 'R_01');
    export_fig([filename '.png']);
    export_fig([filename '.pdf']);
    saveas(gcf, filename, 'fig')
end
%%
X1 = U(pts, pts);
X1 = X1 ./ max(abs(X1(:)));
X1 = X1(inds, :);
h1 = schemaball(X1);
title('Just U')
if SAVE
    filename = fullfile(out_folder2, 'U_00');
    export_fig([filename '.png']);
    export_fig([filename '.pdf']);
    saveas(gcf, filename, 'fig')
end

%% Plot
figure
subplot(131)
plot_resistance_u(m, R, u, VERT)
title({'R', sprintf('elapsed=%g', elapsed_R)})
subplot(132)
plot_resistance_u(m, squareform(Rtilde), u, VERT)
title({'Rtilde', sprintf('elapsed=%g', elapsed_Rtilde)})
subplot(133)

figure
subplot(331) ; title('resistance_{ij} - u_i + u_j') ; axis off
subplot(332) ; plot_resistance_u(m, R, u, VERT) ; title('R') ; colorbar
subplot(333) ; plot_resistance_u(m, squareform(Rtilde), u ,VERT) ; title('Rtilde') ; colorbar
subplot(334) ; title('resistance_{ij}') ;  axis off
subplot(335) ; plot_resistance(m, R, VERT) ; colorbar
subplot(336) ; plot_resistance(m, squareform(Rtilde), VERT); colorbar
subplot(337) ; title('-u_i + u_j') ; axis off
subplot(338) ; m.draw(-u(VERT)+u, 'FaceAlpha', 1, 'EdgeAlpha', 0) ; colorbar
subplot(339) ; m.draw(-u(VERT)+u, 'FaceAlpha', 1, 'EdgeAlpha', 0) ; colorbar


if SAVE
    filename = fullfile(out_folder, [meshname, '_R_vs_Rtilde']);
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end

%%
Rtilde = {R};
names = {{'R', sprintf('elapsed=%g', elapsed_R)}};
i = 2;
progressbar
for eps = 0.2:0.1:1
    rng(SEED)
    tic
    Rtilde{i} = resistance_distance(m, eps, JLFAC, USE_GPU); 
    if ~isempty(gd), wait(gd); end; elapsed = toc;
    names{i} = {sprintf('eps=%g', eps), ...
        sprintf('err=%g', norm(R-squareform(Rtilde{i}), 'fro')) / norm(R, 'fro'), ...
        sprintf('elapsed=%g', elapsed)};
    progressbar((i-1) / length(0.2:0.1:1))
    i = i + 1;
end

figure
for i = 1:numel(names)
    subplot(3, ceil(numel(names) / 3), i)
    X = Rtilde{i};
    if size(X, 2) == 1
        X = squareform(X);
    end
    plot_resistance_u(m, X, u)
    title(names{i})
end

if SAVE
    filename = fullfile(out_folder, [meshname, '_R_vs_Rtilde_varying_eps']);
    print(gcf, filename, '-dpng', '-r300')
    saveas(gcf, filename, 'fig')
end