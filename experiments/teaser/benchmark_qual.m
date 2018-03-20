%% torus_fat_r2
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex')
    
%% pensatore
close all

DEGREE         = 4;
LW             = 15;
MAX_ANIM_T     = 50;
PERCENTAGE     = 0.005;
DATA_FOLDER    = fullfile('..', '..', '..', 'data', 'ashish_nob');
FFIELD_FOLDER  = fullfile('..', 'ioq_benchmark', 'results', 'ashish_nob_02', 'ffields');
STATS_FILE     = fullfile('..', 'ioq_benchmark', 'results', 'ashish_nob_02', 'stats');
BIN_PATH       = fullfile('..', '..', 'vis_vector_field', 'vis_vector_field_bin.exe');
CAMERA_ZOOM    = 1.61675;
FONT_SIZE      = 14;
SAVE_TITLE     = true;
CLOSE          = true;

% name = 'kitten100K'; %Pretty good
% MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);
% CM              = bone; % goes well with mesh color .2.2.4
% MESH_COLOR      = '1,1,1';
% POS_COLOR       = '0,.8,.2';
% NEG_COLOR       = '1,0,0';
% SING_SIZE      = 120;
% % model angle (a quaternion, it's printed to the cmd so you can take it from there)
% ANGLEs          = {'-0.0873874,-0.0961444,-0.00847384,-0.991488',...; % front
%                    '-0.0178393,0.368089,0.00706389,-0.929593',...
%                    '-0.000399498,0.99338,0.00345714,-0.114821'}; % back
% angle_names = {'front','side','back'};               

% name = 'pensatore'; %Pretty good
% MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);
% CM              = bone; % goes well with mesh color .2.2.4
% cmm = pink;  
% MESH_COLOR      = '.72,.68,1';
% POS_COLOR       = '0,.8,.2';
% NEG_COLOR       = '1,0,0';
% SING_SIZE      = 120;
% % model angle (a quaternion, it's printed to the cmd so you can take it from there)
% ANGLEs          = {'0.179615,0.383047,0.0582874,-0.904221',...; % front
%                    '0.112393,-0.903726,-0.280931,-0.302862 ',... % back
%                    '0.13735,-0.488736,-0.169798,-0.844654',...}; % side
%                    '-0.000892228,0.941778,0.246346,-0.228839'}; % side 2
% angle_names = {'front','back','side','side2'};               

% name = 'elephant'; %Pretty good
% MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);
% CM              = bone; % goes well with mesh color .2.2.4
% MESH_COLOR      = '.6,.8,1';
% POS_COLOR       = '0,.8,.2';
% NEG_COLOR       = '1,0,0';
% SING_SIZE      = 80;
% % model angle (a quaternion, it's printed to the cmd so you can take it from there)
% ANGLEs          = {...
%                    '-0.286883,-0.126014,-0.0891663,0.945446',... % side
%                    '-0.020046,-0.958322,-0.276379,0.0695078',... % side2
%                    '-0.134325,-0.778741,-0.194954,0.580958',...  % front
%                    '0.141554,0.474814,0.403799,0.769064',...     % top
%                    };
% angle_names = {'side','side2','front','top'};               
     
name = 'fertility_tri'; %Pretty good
MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);
CM              = bone; % goes well with mesh color .2.2.4
%MESH_COLOR      = '.9,.7,1';
MESH_COLOR      = '1,.9,.7';
POS_COLOR       = '0,.8,.2';
NEG_COLOR       = '1,0,0';
SING_SIZE      = 80;
% model angle (a quaternion, it's printed to the cmd so you can take it from there)
ANGLEs          = {...
                   '0.196401,0.0871462,-0.171236,0.961515',... % side
                   '0.180944,-0.963498,-0.190592,-0.0510479',... % side2
                   '0.0571014,-0.708373,0.0451759,-0.702073',...  % back
                   };
angle_names = {'side','side2','back'};               

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
       title(sprintf('%s, $E = %.2f, |S| = %d$', alg_title{i}, m.miq_energy, m.n_vert_sing), 'FontSize', FONT_SIZE);
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
