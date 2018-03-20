set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex')

FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = true;
SAVE = false;
LOAD = false;
OUT_FOLDER = 'results';
mkdir(OUT_FOLDER);
DATA_FOLDER = '../../../data/high_genus/';
AX_FONT_SIZE = 10;
FONT_SIZE = 13;

files = get_filepaths(DATA_FOLDER, '.off');

Ea = [];
Eb_rnd = [];
Eab_rnd = [];
Eb_cvp = [];
Eab_cvp = [];

Z_MAX = 20;

progressbar
beta_rnd_results = {};
beta_cvp_results = {};
alpha_ioq_results = {};

T_rnd = [];
T_cvp = [];

for i = 1:numel(files)
    fp = files{i};
    [~, meshname, ~] = fileparts(fp);
    print_header(meshname);
    m = Mesh(fp);
    V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;
    alpha_g = get_gaussian_curvature(m);
    beta_g = wrapToPi(generator_angle_defects(m));
    
    rng(SEED)
    [alpha_ioq, beta_rnd, elapsed_ioq] = IOQ_highgenus_gpu(...
        V, F, ...
        'UseGPU', USE_GPU, ...
        'Iterations', 2000);
    
    d0 = get_exterior_derivatives(m);
    L = d0'*d0;
    H = -m.H;
    %H = [H(:, 1), -H(:, 2)];
    %H = [-H(:, 1), H(:, 2)];
    F = factorize(L);
    A2 = (H'*d0) / F;
    y0 = (2/pi)*(A2*alpha_g - beta_g);
    B = H - d0 * (F \ (d0' * H));
    G = B * inverse(H'*B) * (pi/2);
    
    calc_a = @(alpha) inverse(L) * ( (pi/2)*alpha - alpha_g );
    calc_b = @(alpha, beta) inverse(H'*B)*( (pi/2)*beta - beta_g - H'*d0*calc_a(alpha));
    
    a_ioq = calc_a(alpha_ioq);
    b_rnd = calc_b(alpha_ioq, beta_rnd);
    
    m = size(G, 2);
    target = B * inverse(H'*B) * ( beta_g + H'*d0*a_ioq );
    tic
    [beta_cvp,~,~,~,~] = SEA_det(m, Z_MAX, target, G, false);
    T_cvp(end+1) = toc;
    
    b_cvp = calc_b(alpha_ioq, beta_cvp);
    
    Ea(end+1)      = norm(d0*a_ioq)^2;
    Eb_rnd(end+1)  = norm(B*b_rnd)^2;
    Eab_rnd(end+1) = Ea(end) + Eb_rnd(end);
    Eb_cvp(end+1)  = norm(B*b_cvp)^2 ;
    Eab_cvp(end+1) = Ea(end) + Eb_cvp(end);
    
    beta_rnd_results{end+1} = beta_rnd;
    beta_cvp_results{end+1} = beta_cvp;
    alpha_ioq_results{end+1} = alpha_ioq;
    
    progressbar(i / numel(files))
end

genus = [];
for i = 1:numel(files)
    fp = files{i};
    m = Mesh(fp);
    genus(end+1) = m.genus;
end

%% Create ffields
for i = 1:numel(files)
    fp = files{i};
    [~, meshname, ~] = fileparts(fp);
    m = Mesh(fp);
    
    alpha_ioq = alpha_ioq_results{i};
    beta_rnd = beta_rnd_results{i};
    beta_cvp = beta_cvp_results{i};
    k_rnd = [alpha_ioq; beta_rnd];
    k_cvp = [alpha_ioq; beta_cvp];
    
    rnd_results{i} = TCODS(m, ...
        'k', k_rnd, ...
        'f0', FACE0, ...
        'theta0', THETA0, ...
        'degree', DEGREE, ...
        'CreateFField', true, ...
        'Duplicate', true, ...
        'gConstraintVec', GVEC);
    filename = fullfile(OUT_FOLDER, [meshname '_rnd.ffield']);
    rnd_results{i}.saveTM(filename);
    
    cvp_results{i} = TCODS(m, ...
        'k', k_cvp, ...
        'f0', FACE0, ...
        'theta0', THETA0, ...
        'degree', DEGREE, ...
        'CreateFField', true, ...
        'Duplicate', true, ...
        'gConstraintVec', GVEC);
    filename = fullfile(OUT_FOLDER, [meshname '_cvp.ffield']);
    cvp_results{i}.saveTM(filename);
end

%% Plot

col = linspecer(3); cc = col(2,:); col(2,:) = col(3,:); col(3,:) = cc;

% sort
[genus2, is] = sort(genus);
Ea2 = Ea(is);
Eb_rnd2 = Eb_rnd(is);
Eab_rnd2 = Eab_rnd(is);
Eb_cvp2 = Eb_cvp(is);
Eab_cvp2 = Eab_cvp(is);

alpha_ioq_results2 = alpha_ioq_results(is);
beta_rnd_results2 = beta_rnd_results(is);
beta_cvp_results2 = beta_cvp_results(is);

close all
names = cellfun(@(x) strrep(get_name(x), '_', '\_'), files, 'UniformOutput', false);
names2 = names(is);
fh = figure;
plot(1:numel(files), Eb_rnd2, 'x', ...
     1:numel(files), Eb_cvp2, 'o');
legend('Round', 'CVP', 'Location', 'northwest')
%plt = Plot(fh, true);
%plt.ShowBox = 'off';
%plt.Markers = {'d', '.', 'x'};
%plt.Legend = {'Ea', 'Eb rnd', 'Eb cvp'};

%plt.LineStyle = {':', ':', ':'};
%plt.Markers = {'d', '.', 'x'};

%%
fh = figure;
c = categorical(names2);
c = reordercats(c, names2);
%h = barh(c,[Eb_rnd2, Eb_cvp2],'edgecolor','none'); 
h = barh(c,[Eb_rnd2', Eb_cvp2'],'edgecolor','none'); 
set(gca, 'TickLabelInterpreter', 'latex');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
lh = legend('Round', 'CVP');
lh.FontSize = FONT_SIZE;
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
ax.FontSize = AX_FONT_SIZE;
plt = Plot(fh, true);
plt.BoxDim = [6, 3];
plt.ShowBox = 'off';
plt.XGrid='off';

export_fig('cvp_benchmark.pdf')

%% Timing
plt = Plot(T_cvp(is));
