%% Compare IOQ conn with IOQ cot

%% Setup
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex')
    
clear all; close all;

name = 'decimated-max';
%row  = 9; 

DEGREE         = 4;
LW             = 15;
MAX_ANIM_T     = 50;
PERCENTAGE     = 0.005;
%DATA_FOLDER   = fullfile('..', '..', '..', 'data', 'ashish_nob');
DATA_FOLDER    = fullfile('..', '..', '..', 'data');
OUT_FOLDER     = fullfile('results', [name, '-vf']);
mkdir(OUT_FOLDER);
%FFIELD_FOLDER = fullfile('..', 'ioq_benchmark', 'results', 'ashish_nob_02', 'ffields');
MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);
BIN_PATH       = fullfile('..', '..', 'vis_vector_field', 'vis_vector_field_bin.exe');
CM              = bone; % goes well with mesh color .2.2.4
% Various colors (make sure there are no spaces)
cmm = pink;  
MESH_COLOR      = '1,1,1';
POS_COLOR       = '0,.8,.2';
NEG_COLOR       = '1,0,0';
SING_SIZE       = 120;
% model angle (a quaternion, it's printed to the cmd so you can take it from there)
ANGLEs          = {'-0.0873874,-0.0961444,-0.00847384,-0.991488',...; % front
                   '-0.0178393,0.368089,0.00706389,-0.929593',...
                   '-0.000399498,0.99338,0.00345714,-0.114821'}; % back
angle_names     = {'front','side','back'};               
CAMERA_ZOOM     = 1.61675;
FONT_SIZE       = 14;
SAVE_TITLE      = true;
CLOSE           = true;
SEED            = 112;
USE_GPU         = true;
INV_METHOD      = 'GPUInv';
FACE0 = 1; THETA0 = 0; GVEC = [1,0,0];

algs = {'IOQ conn','IOQ cot'};
%alg_title = {'MIQ','GO','IOQ','IOQe $\epsilon=.5$'};

%algs = {'MIQ'};
%cols = [2, 3, 4, 6];

vis_props = {'CM', CM, ...
        'degree', DEGREE, ...
        'lw', LW, ...
        'max_anim_t', MAX_ANIM_T, ...
        'percentage', PERCENTAGE, ...
        'mesh_color', MESH_COLOR, ...
        'pos_color', POS_COLOR, ...
        'neg_color', NEG_COLOR, ...
        'sing_size', SING_SIZE, ...
        'angle', ANGLE, ...
        'camera_zoom', CAMERA_ZOOM, ...
        'Close', CLOSE, ...            
        'BinaryPath', BIN_PATH};

%% Load mesh
m = Mesh(MESH_FP);
nv = m.nV; nf = m.nF; ne = m.nE;
V = m.V; F = m.F;

%% Run IOQ conn

rng(SEED)
[alpha1, beta1, elapsed_ioq1] = IOQ_highgenus_gpu(...
    V, F, ...
    'UseGPU', USE_GPU, ...
    'Iterations', 2000, ...
    'Laplacian', 'conn', ...
    'InvMethod', INV_METHOD);

%
k1 = [alpha1; beta1];
res_ioq_conn = TCODS(m, ...
    'k', k1, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);


%% run IOQ with cot lap
rng(SEED)
[alpha2, beta2, elapsed_ioq2] = IOQ_highgenus_gpu(...
    V, F, ...
    'UseGPU', USE_GPU, ...
    'Iterations', 2000, ...
    'Laplacian', 'cot', ...
    'InvMethod', INV_METHOD);

%
k2 = [alpha2; beta2];
res_ioq_cot = TCODS(m, ...
    'k', k2, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);


%% run MIQ
[res_miq, elapsed_miq] = ...
    nrosy_mex(MESH_FP, FACE0, GVEC, DEGREE);


%% Visualize

title_ioq_conn = {'', sprintf('IOQ $E = %g, |S| = %d\\hskip 4em$',      res_ioq_conn.miq_energy, res_ioq_conn.n_vert_sing)};
title_ioq_cot = {'', sprintf('IOQc $E = %g, |S| = %d$\\hskip 4em', res_ioq_cot.miq_energy,  res_ioq_cot.n_vert_sing)};
title_miq = {'', sprintf('MIQ $E = %g, |S| = %d$',      res_miq.miq_energy,      res_ioq_conn.n_vert_sing)};

ANGLE = ANGLEs{3};
names = {};

R = res_ioq_conn.ffield_vectors(1:nf, :);
outname = fullfile(OUT_FOLDER, 'ioq_conn.png');
vis_vector_field(MESH_FP, ...
        R, ...
        alpha1, ...
        'out', outname, ...
        vis_props{:});
outname_title = fullfile(OUT_FOLDER, 'ioq_conn_title.png');
create_title(outname_title, title_ioq_conn, FONT_SIZE);
out_pan = fullfile(OUT_FOLDER, 'ioq_conn2.png');
panorama({outname, outname_title}, out_pan, 'vertical');
names{end+1} = out_pan;

R = res_ioq_cot.ffield_vectors(1:nf, :);
outname = fullfile(OUT_FOLDER, 'ioq_cot.png');
vis_vector_field(MESH_FP, ...
        R, ...
        alpha2, ...
        'out', outname, ...
        vis_props{:});
outname_title = fullfile(OUT_FOLDER, 'ioq_cot_title.png');
create_title(outname_title, title_ioq_cot, FONT_SIZE);
out_pan = fullfile(OUT_FOLDER, 'ioq_cot2.png');
panorama({outname, outname_title}, out_pan, 'vertical');
names{end+1} = out_pan;

R = res_miq.ffield_vectors(1:nf, :);
alpha3 = zeros(nv, 1);
inds = res_miq.vert_sing(:, 1); alpha3(inds) = DEGREE*res_miq.vert_sing(:,2);
outname = fullfile(OUT_FOLDER, 'miq.png');
vis_vector_field(MESH_FP, ...
        R, ...
        alpha3, ...
        'out', outname, ...
        vis_props{:});
outname_title = fullfile(OUT_FOLDER, 'miq_title.png');
create_title(outname_title, title_miq, FONT_SIZE);
out_pan = fullfile(OUT_FOLDER, 'miq2.png');
panorama({outname, outname_title}, out_pan, 'vertical');
names{end+1} = out_pan;

out_pan = fullfile(OUT_FOLDER, [name, '_pan.png']);
panorama(names, out_pan);


% for i = 1:length(algs)
%     outnames = {};
% 
%     for j = 1:length(ANGLEs)
%         ANGLE = ANGLEs{j};
%         aname = angle_names{j};
%         
%         alg = algs{i};
%         name1 = [name '_' alg];
% 
%         outname = [name1 '_' aname '.png']
%         outnames{end+1} = outname;
%         
%         ffield_fp1 = fullfile(FFIELD_FOLDER, [name1, '.ffield']);
%         m = Mesh(ffield_fp1); nf = m.nF; nv = m.nV;
% 
%         R = m.ffield_vectors(1:nf, :);
%         S = zeros(nv, 1);
%         S(m.vert_sing(:, 1)) = m.vert_sing(:, 2);
% 
%         vis_vector_field(MESH_FP, ...
%                 R, ...
%                 S, ...
%                 'CM', CM, ...
%                 'degree', DEGREE, ...
%                 'lw', LW, ...
%                 'max_anim_t', MAX_ANIM_T, ...
%                 'percentage', PERCENTAGE, ...
%                 'mesh_color', MESH_COLOR, ...
%                 'pos_color', POS_COLOR, ...
%                 'neg_color', NEG_COLOR, ...
%                 'out', outname, ...
%                 'sing_size', SING_SIZE, ...
%                 'angle', ANGLE, ...
%                 'camera_zoom', CAMERA_ZOOM, ...
%                 'Close', CLOSE, ...            
%                 'BinaryPath', BIN_PATH);
%     end
%     if SAVE_TITLE
%         elapsed = T(row, cols(i));
%         sprintf('%s $E = %.2f, ns = %d, T = %.2f$', name1, m.miq_energy, m.n_vert_sing, elapsed)
%        figure
% %        title(sprintf('$E = %.2f, |S| = %d, T = %.2f$', m.miq_energy, m.n_vert_sing, elapsed));
%        title(sprintf('%s, $E = %.2f, |S| = %d$', alg_title{i}, m.miq_energy, m.n_vert_sing), 'FontSize', FONT_SIZE);
%         axis off
%         set(gcf, 'WindowStyle', 'docked')
%         set(gcf,'color','w');
% 
%        export_fig([name1 '_title.pdf']);
%     end
%     
%     out_pan = [name1 '.png'];
%     pannames{end+1} = out_pan;
%     cmd = ['"C:/Program Files (x86)/IrfanView/i_view32.exe" /panorama=(1'];
%     for j=1:length(ANGLEs)
%         cmd = [cmd ',' outnames{j}];
%     end
%     cmd = [cmd ') /convert=' out_pan];
%     system(cmd)    
% end
% 
% out_pan = [name '.png'];
% cmd = ['"C:/Program Files (x86)/IrfanView/i_view32.exe" /panorama=(2'];
% for i=1:length(algs)
%     cmd = [cmd ',' pannames{i}];
% end
% cmd = [cmd ') /convert=' out_pan];
% system(cmd)

%%
set(groot, 'defaultAxesTickLabelInterpreter','none'); 
set(groot, 'defaultLegendInterpreter','none');
set(0,'defaulttextinterpreter','none')
