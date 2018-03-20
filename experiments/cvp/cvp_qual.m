%% torus_fat_r2
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex')
    
%% pensatore
close all

DEGREE         = 4;
LW             = 15;
MAX_ANIM_T     = 50;
DATA_FOLDER    = fullfile('..', '..', '..', 'data', 'high_genus');
FFIELD_FOLDER  = fullfile('..', 'cvp', 'results');
BIN_PATH       = fullfile('..', '..', 'vis_vector_field', 'vis_vector_field_bin.exe');
CAMERA_ZOOM    = 1.61675;
FONT_SIZE      = 14;
SAVE_TITLE     = false;
CLOSE          = true;

name = 'ball_r'; %Pretty good
MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);
CM             = bone; % goes well with mesh color .2.2.4
MESH_COLOR     = '1,1,1';
POS_COLOR      = '0,.8,.2';
NEG_COLOR      = '1,0,0';
SING_SIZE      = 120;
PERCENTAGE     = 0.05;
% model angle (a quaternion, it's printed to the cmd so you can take it from there)
ANGLEs          = {'-0.283409,-0.363991,-0.117299,-0.879449',...; % front
                   }; % back
angle_names = {'front'};               

% name = 'torus_s0'; %Pretty good
% MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);
% CM              = bone; % goes well with mesh color .2.2.4
% MESH_COLOR      = '1,1,1';
% POS_COLOR       = '0,.8,.2';
% NEG_COLOR       = '1,0,0';
% SING_SIZE      = 120;
% PERCENTAGE     = 0.1;
% % model angle (a quaternion, it's printed to the cmd so you can take it from there)
% ANGLEs          = {'0.619045,-0.333941,0.32774,-0.630757',...; % front
%                    }; % back
% angle_names = {'front'};               

% name = 'dancing_children_r'; %Pretty good
% MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);
% CM              = bone; % goes well with mesh color .2.2.4
% MESH_COLOR      = '1,1,1';
% POS_COLOR       = '0,.8,.2';
% NEG_COLOR       = '1,0,0';
% SING_SIZE      = 30;
% PERCENTAGE     = 0.1;
% % model angle (a quaternion, it's printed to the cmd so you can take it from there)
% ANGLEs          = {'0.0777115,0.0876955,0.00686262,0.993088',...; % front
%                    '-0.00391143,-0.988324,-0.150121,-0.0257525 ',... % back
%                    }; % back
% angle_names = {'front','back'};               

% name = 'bob_tri'; %Pretty good
% MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);
% CM              = cool; % goes well with mesh color .2.2.4
% MESH_COLOR      = '.2,.2,.4';
% POS_COLOR       = '0,.8,.2';
% NEG_COLOR       = '1,0,0';
% SING_SIZE      = 30;
% PERCENTAGE     = 0.05;
% % model angle (a quaternion, it's printed to the cmd so you can take it from there)
% ANGLEs          = {'-0.423573,-0.310589,-0.157312,-0.836285',...; % front
%                    }; % back
% angle_names = {'front'};     

% name = 'eight'; %Pretty good
% MESH_FP        = fullfile(DATA_FOLDER, [name '.off']);
% CM              = cool; % goes well with mesh color .2.2.4
% MESH_COLOR      = '.2,.2,.4';
% POS_COLOR       = '0,.8,.2';
% NEG_COLOR       = '1,0,0';
% SING_SIZE      = 30;
% PERCENTAGE     = 0.05;
% % model angle (a quaternion, it's printed to the cmd so you can take it from there)
% ANGLEs          = {'-0.423573,-0.310589,-0.157312,-0.836285',...; % front
%                    }; % back
% angle_names = {'front'};     



algs = {'rnd','cvp'};
alg_title = algs;
alg_fname = algs;

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
    %cmd = ['"C:/Program Files/IrfanView/i_view64.exe" /panorama=(1'];
    for j=1:length(ANGLEs)
        cmd = [cmd ',' outnames{j}];
    end
    cmd = [cmd ') /resize=(30p,30p) /resample /convert=' out_pan];
    system(cmd)    
end

out_pan = [name '.jpg'];
cmd = ['"C:/Program Files (x86)/IrfanView/i_view32.exe" /panorama=(1'];
for i=1:length(algs)
    cmd = [cmd ',' pannames{i}];
end
cmd = [cmd ') /convert=' out_pan];
system(cmd)

%%
set(groot, 'defaultAxesTickLabelInterpreter','none'); 
set(groot, 'defaultLegendInterpreter','none');
set(0,'defaulttextinterpreter','none')
