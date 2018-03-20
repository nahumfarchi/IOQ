%% torus_fat_r2
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex')
    
%% pensatore
close all

%name = 'bimba100K'; %Pretty good
name = 'bunny';
%row  = 9; 
%outdir = ['results/' name '/'];
outdir = './';
if ~exist(outdir) 
    mkdir(outdir);
end

DEGREE         = 4;
LW             = 15;
MAX_ANIM_T     = 50;
PERCENTAGE     = 0.005;
%DATA_FOLDER   = fullfile('..', '..', '..', 'data', 'ashish_nob');
%DATA_FOLDER    = fullfile('..', '..', '..', 'data', 'ashish_nob');
%FFIELD_FOLDER = fullfile('..', 'ioq_benchmark', 'results', 'ashish_nob_02', 'ffields');
%FFIELD_FOLDER  = fullfile('..', 'ioq_benchmark', 'results', 'ashish_nob_02', 'ffields');
DATA_FOLDER    = './';
FFIELD_FOLDER  = './';
STATS_FILE     = fullfile('..', 'ioq_benchmark', 'results', 'ashish_nob_02', 'stats');
MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);
BIN_PATH       = fullfile('..', '..', 'vis_vector_field', 'vis_vector_field_bin.exe');
CM             = cool; % goes well with .2.2.4
%CM              = winter;
%CM              = bone; % goes well with mesh color .2.2.4
%CM              = pink; % goes well with mesh color .2.2.4
%CM              = copper; % goes well with mesh color '.86,.79,.64'
%CM             = zeros(40, 3);
% Various colors (make sure there are no spaces)
cmm = pink;  
%MESH_COLOR     = '.86,.79,.64';%'.72,.68,1';
MESH_COLOR     = '.2,.2,.4';
%POS_COLOR      = '0.2157,0.4941,0.7216';
POS_COLOR       = '0,.8,.2';
%NEG_COLOR      = '0.8941,0.1020,0.1098';
NEG_COLOR       = '1,0,0';
SING_SIZE      = 120;
% model angle (a quaternion, it's printed to the cmd so you can take it from there)
% ANGLEs          = {'0.121267,0.311376,-0.0400987,-0.941664',...; % Front
%                    '-0.00314773,0.989292,0.144311,-0.0215794',... % Back
%                    '-0.325132,0.507198,0.214503,-0.768784'}; % top
ANGLEs = {[], [], []};
angle_names = {'front','back','top'};               
%CAMERA_ZOOM    = 1.61675;
CAMERA_ZOOM = [];
FONT_SIZE      = 14;
SAVE_TITLE     = true;
CLOSE          = false;


s = load(STATS_FILE); s = s.stats;
T = cellfun(@get_time, s, 'UniformOutput', true);
row = find(cellfun(@(x) strcmp(x, name), s(:, 1)));

%algs = {'MIQ','IOQ_conn_1a','JL_IOQ_eps0.5','GO'};
%algs = {'MIQ','GO','IOQ_conn_1a', 'JL_IOQ_eps0.5'};
%algs = {'MIQ','JL_IOQ_eps0.5'};
%alg_title = {'MIQ','IOQe $\epsilon=.5$'};
algs = {'01'};
alg_title = {'IOQe'};
%algs = {'MIQ'};
cols = [2, 6];

pannames = {};

for i = 1:length(algs)
    outnames = {};
    
    for j = 1:length(ANGLEs)
        ANGLE = ANGLEs{j};
        aname = angle_names{j};
        
        alg = algs{i};
        name1 = [name '_' alg];

        outname = [name1 '_' aname '.png']
        outnames{end+1} = outname;
        
        ffield_fp1 = fullfile(FFIELD_FOLDER, [name1, '.ffield']);
        m = Mesh(ffield_fp1); nf = m.nF; nv = m.nV;

        R = m.ffield_vectors(1:nf, :);
        S = zeros(nv, 1);
        S(m.vert_sing(:, 1)) = m.vert_sing(:, 2);

        MESH_FP = [name '_' alg '.off'];
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
                'out', outname, ...
                'sing_size', SING_SIZE, ...
                'angle', ANGLE, ...
                'camera_zoom', CAMERA_ZOOM, ...
                'Close', CLOSE, ...            
                'BinaryPath', BIN_PATH);
    end
    if SAVE_TITLE
        elapsed = T(row, cols(i));
        sprintf('%s $E = %.2f, ns = %d, T = %.2f$', name1, m.miq_energy, m.n_vert_sing, elapsed)
       figure
%        title(sprintf('$E = %.2f, |S| = %d, T = %.2f$', m.miq_energy, m.n_vert_sing, elapsed));
       title(sprintf('%s, $E = %.2f, |S| = %d$', alg_title{i}, m.miq_energy, m.n_vert_sing), 'FontSize', FONT_SIZE);
        axis off
        set(gcf, 'WindowStyle', 'docked')
        set(gcf,'color','w');

       export_fig([outdir name1 '_title.pdf']);
    end
    
    out_pan = [name1 '.png'];
    pannames{end+1} = out_pan;
    cmd = ['"C:/Program Files (x86)/IrfanView/i_view32.exe" /panorama=(1'];
    for j=1:length(ANGLEs)
        cmd = [cmd ',' outnames{j}];
    end
    cmd = [cmd ') /convert=' out_pan];
    system(cmd)
end

out_pan = [name '.png'];
cmd = ['"C:/Program Files (x86)/IrfanView/i_view32.exe" /panorama=(2'];
for i=1:length(algs)
    cmd = [cmd ',' pannames{i}];
end
cmd = [cmd ') /convert=' out_pan];
system(cmd)

%%
set(groot, 'defaultAxesTickLabelInterpreter','none'); 
set(groot, 'defaultLegendInterpreter','none');
set(0,'defaulttextinterpreter','none')
