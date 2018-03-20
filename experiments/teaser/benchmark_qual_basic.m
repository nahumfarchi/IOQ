set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex')
    
close all
clear all



% name = 'rounded_cube_keenan';
% dataset = 'basic';
% dataset_exp = dataset;
% CM              = cool; % goes well with mesh color .2.2.4
% MESH_COLOR      = '.2,.2,.4';
% POS_COLOR       = '0,.8,.2';
% NEG_COLOR       = '1,0,0';
% SING_SIZE      = 180;
% % model angle (a quaternion, it's printed to the cmd so you can take it from there)
% ANGLEs          = {...
%                    '-0.102709,0.855176,0.298207,-0.411336'
%                    };
% CAMERA_ZOOM    = 1.1;               
% angle_names = {'front'};               
%      
name = 'torus_fat_r2';
dataset = 'genus1_small';
dataset_exp = [dataset '_02'];
CM              = bone; % goes well with mesh color .2.2.4
MESH_COLOR      = '.2,0,.1';
POS_COLOR       = '0,.8,.2';
NEG_COLOR       = '1,0,0';
SING_SIZE      = 180;
% model angle (a quaternion, it's printed to the cmd so you can take it from there)
ANGLEs          = {...
                   '-0.690327,-0.196511,-0.219915,-0.660658'
                   };
CAMERA_ZOOM    = 1.61675;
angle_names = {'front'};       

DEGREE         = 4;
LW             = 15;
MAX_ANIM_T     = 50;
PERCENTAGE     = 0.005;
DATA_FOLDER    = fullfile('..', '..', '..', 'data', dataset);
FFIELD_FOLDER  = fullfile('..', 'ioq_benchmark', 'results', dataset_exp, 'ffields');
STATS_FILE     = fullfile('..', 'ioq_benchmark', 'results', dataset_exp, 'stats');
BIN_PATH       = fullfile('..', '..', 'vis_vector_field', 'vis_vector_field_bin.exe');
FONT_SIZE      = 14;
SAVE_TITLE     = true;
CLOSE          = true;
MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);

s = load(STATS_FILE); s = s.stats;
T = cellfun(@get_time, s, 'UniformOutput', true);
row = find(cellfun(@(x) strcmp(x, name), s(:, 1)));

algs = {'MIQ','GO','IOQ_conn_1a','JL_IOQ_eps0.5'};
alg_title = {'MIQ','GO','IOQ','IOQe $\epsilon=.5$'};
alg_fname = {'MIQ','GO','IOQ','IOQe0_5'};

cols = [2, 3, 4, 6];

pannames = {};

for i = 1:length(algs)
    outnames = {};

    for j = 1:length(ANGLEs)
        ANGLE = ANGLEs{j};
        aname = angle_names{j};
        
        alg = algs{i};
        name1 = [name '_' algs{i}];

        outname = [name1 '_' aname '.png']
        outnames{end+1} = outname;
        
        ffield_fp1 = fullfile(FFIELD_FOLDER, [name '_' alg '.ffield']);
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
%        title(sprintf('%s, $E = %.2f, |S| = %d$\n$T = %.2f$ s', ...
%            alg_title{i}, m.miq_energy, m.n_vert_sing,elapsed), 'FontSize', FONT_SIZE,'horizontalAlignment', 'left');
       title(sprintf('%s, $E = %.2f, |S| = %d$', ...
           alg_title{i}, m.miq_energy, m.n_vert_sing), 'FontSize', FONT_SIZE);
         axis off
        set(gcf, 'WindowStyle', 'docked')
        set(gcf,'color','w');

       export_fig([name '_' alg_fname{i} '_title.pdf']);
    end
    
    out_pan = [name '_' alg_fname{i} '.jpg'];
    pannames{end+1} = out_pan;
    cmd = ['"C:/Program Files (x86)/IrfanView/i_view32.exe" /panorama=(1'];
    for j=1:length(ANGLEs)
        cmd = [cmd ',' outnames{j}];
    end
    cmd = [cmd ') /resize=(30p,30p) /resample /convert=' out_pan];
    system(cmd)    
end

out_pan = [name '.jpg'];
pandir = '2'; % vertical
if length(ANGLEs) == 1
    pandir = '1'; %horizontal
end
cmd = ['"C:/Program Files (x86)/IrfanView/i_view32.exe" /panorama=(' num2str(pandir)];
for i=1:length(algs)
    cmd = [cmd ',' pannames{i}];
end
cmd = [cmd ') /convert=' out_pan];
system(cmd)

%%
set(groot, 'defaultAxesTickLabelInterpreter','none'); 
set(groot, 'defaultLegendInterpreter','none');
set(0,'defaulttextinterpreter','none')
