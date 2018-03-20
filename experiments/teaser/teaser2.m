%% torus_fat_r2
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex')

close all

name = 'torus_fat_r2';
row  = 9; 

DEGREE         = 4;
LW             = 2;
MAX_ANIM_T     = 50;
PERCENTAGE     = 0.01;
%DATA_FOLDER   = fullfile('..', '..', '..', 'data', 'ashish_nob');
DATA_FOLDER    = fullfile('..', '..', '..', 'data', 'genus1_small');
%FFIELD_FOLDER = fullfile('..', 'ioq_benchmark', 'results', 'ashish_nob_02', 'ffields');
FFIELD_FOLDER  = fullfile('..', 'ioq_benchmark', 'results', 'genus1_small_02', 'ffields');
STATS_FILE     = fullfile('..', 'ioq_benchmark', 'results', 'genus1_small_02', 'stats');
MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);
BIN_PATH       = fullfile('..', '..', 'vis_vector_field', 'vis_vector_field_bin.exe');
CM             = linspecer(40);
% Various colors (make sure there are no spaces)
MESH_COLOR     = '0,0,0';
POS_COLOR      = '0.2157,0.4941,0.7216';
NEG_COLOR      = '0.8941,0.1020,0.1098';
SING_SIZE      = 120;
% model angle (a quaternion, it's printed to the cmd so you can take it from there)
%ANGLE          = '0.230011,0.711119,0.607353,0.269308';
ANGLE          = '-0.101072,0.717077,-0.681374,0.106369';
CAMERA_ZOOM    = 1.88565;
FONT_SIZE      = 14;
SAVE_TITLE     = true;
CLOSE          = true;

s = load(STATS_FILE); s = s.stats;
T = cellfun(@get_time, s, 'UniformOutput', true);

%algs = {'MIQ','IOQ_conn_1a','JL_IOQ_eps0.5','GO'};
algs = {'MIQ','IOQ_conn_1a', 'JL_IOQ_eps0.5'};
%algs = {'MIQ'};
cols = [2, 4, 6];

print_header('press any key to save to png')

for i = 1:length(algs)
    
    alg = algs{i};
    name1 = [name '_' alg];

    ffield_fp1 = fullfile(FFIELD_FOLDER, [name1, '.ffield']);
    m = Mesh(ffield_fp1); nf = m.nF; nv = m.nV;

    R = m.ffield_vectors(1:nf, :);
    S = zeros(nv, 1);
    S(m.vert_sing(:, 1)) = m.vert_sing(:, 2);

    vis_vector_field(MESH_FP, ...
            R, ...
            S, ...
            'CM', CM, ...
            'degree', DEGREE, ...
            'lw', LW, ...
            'max_anim_t', MAX_ANIM_T, ...
            'percentage', PERCENTAGE, ...
            'mesh_color', MESH_COLOR, ...
            'pos_color', POS_COLOR, ...
            'neg_color', NEG_COLOR, ...
            'out', [name1 '.png'], ...
            'sing_size', SING_SIZE, ...
            'angle', ANGLE, ...
            'camera_zoom', CAMERA_ZOOM, ...
            'Close', CLOSE, ...
            'BinaryPath', BIN_PATH);
    
    if SAVE_TITLE
        elapsed = T(row, cols(i));
        figure
        %title(sprintf('$E = %.2f, ns = %d, T = %.3f$', m.miq_energy, m.n_vert_sing, elapsed));
        title(sprintf('$E = %.2f, ns = %d$', m.miq_energy, m.n_vert_sing), 'FontSize', FONT_SIZE);
        axis off
        set(gcf, 'WindowStyle', 'docked')
        set(gcf,'color','w');

        export_fig([name1 '_title.pdf']);
    end
    
end
    
%% pensatore
close all

name = 'pensatore';
%row  = 9; 

DEGREE         = 4;
LW             = 2;
MAX_ANIM_T     = 50;
PERCENTAGE     = 0.01;
%DATA_FOLDER   = fullfile('..', '..', '..', 'data', 'ashish_nob');
DATA_FOLDER    = fullfile('..', '..', '..', 'data', 'ashish_nob');
%FFIELD_FOLDER = fullfile('..', 'ioq_benchmark', 'results', 'ashish_nob_02', 'ffields');
FFIELD_FOLDER  = fullfile('..', 'ioq_benchmark', 'results', 'ashish_nob_02', 'ffields');
STATS_FILE     = fullfile('..', 'ioq_benchmark', 'results', 'ashish_nob_02', 'stats');
MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);
BIN_PATH       = fullfile('..', '..', 'vis_vector_field', 'vis_vector_field_bin.exe');
%CM             = linspecer(40);
CM             = zeros(40, 3);
% Various colors (make sure there are no spaces)
%MESH_COLOR     = '0,0,0';
MESH_COLOR     = '1,1,1';
POS_COLOR      = '0.2157,0.4941,0.7216';
NEG_COLOR      = '0.8941,0.1020,0.1098';
SING_SIZE      = 120;
% model angle (a quaternion, it's printed to the cmd so you can take it from there)
ANGLE          = '0.121267,0.311376,-0.0400987,-0.941664';
CAMERA_ZOOM    = 1.61675;
FONT_SIZE      = 14;
SAVE_TITLE     = false;
CLOSE          = true;

s = load(STATS_FILE); s = s.stats;
T = cellfun(@get_time, s, 'UniformOutput', true);
row = find(cellfun(@(x) strcmp(x, name), s(:, 1)));

%algs = {'MIQ','IOQ_conn_1a','JL_IOQ_eps0.5','GO'};
algs = {'MIQ','IOQ_conn_1a', 'JL_IOQ_eps0.5'};
%algs = {'MIQ'};
cols = [2, 4, 6];

print_header('press any key to save to png')

for i = 1:length(algs)
    
    alg = algs{i};
    name1 = [name '_' alg];

    ffield_fp1 = fullfile(FFIELD_FOLDER, [name1, '.ffield']);
    m = Mesh(ffield_fp1); nf = m.nF; nv = m.nV;

    R = m.ffield_vectors(1:nf, :);
    S = zeros(nv, 1);
    S(m.vert_sing(:, 1)) = m.vert_sing(:, 2);

    vis_vector_field(MESH_FP, ...
            R, ...
            S, ...
            'CM', CM, ...
            'degree', DEGREE, ...
            'lw', LW, ...
            'max_anim_t', MAX_ANIM_T, ...
            'percentage', PERCENTAGE, ...
            'mesh_color', MESH_COLOR, ...
            'pos_color', POS_COLOR, ...
            'neg_color', NEG_COLOR, ...
            'out', [name1 '.png'], ...
            'sing_size', SING_SIZE, ...
            'angle', ANGLE, ...
            'camera_zoom', CAMERA_ZOOM, ...
            'Close', CLOSE, ...
            'BinaryPath', BIN_PATH);
        
    if SAVE_TITLE
        elapsed = T(row, cols(i));
        figure
        %title(sprintf('$E = %.2f, ns = %d, T = %.2f$', m.miq_energy, m.n_vert_sing, elapsed));
        title(sprintf('$E = %.2f, ns = %d$', m.miq_energy, m.n_vert_sing), 'FontSize', FONT_SIZE);
        axis off
        set(gcf, 'WindowStyle', 'docked')
        set(gcf,'color','w');
    
        export_fig([name1 '_title.pdf']);
    end
    
end

%%
set(groot, 'defaultAxesTickLabelInterpreter','none'); 
set(groot, 'defaultLegendInterpreter','none');
set(0,'defaulttextinterpreter','none')
