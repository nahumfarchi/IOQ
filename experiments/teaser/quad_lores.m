clear all
close all 

COMPUTE = true;
VIS = false;
% meshname      = 'bimba100K';
% ffield_folder = '../ioq_benchmark/results/ashish_nob_02/ffields';

%meshname      = 'elephant';
%ffield_folder = '../ioq_benchmark/results/ashish_nob_02/ffields';

%meshname      = 'elephant_r';
%ffield_folder = '../ioq_benchmark/results/genus1_small_02/ffields/';

meshname      = 'cat';
ffield_folder = '../ioq_benchmark/results/basic/ffields/';

% meshname      = 'ramses';
% ffield_folder = '../ioq_benchmark/results/ashish_nob_02/ffields/';

%meshname      = 'kitten100k';
%ffield_folder = '../ioq_benchmark/results/ashish_nob_02/ffields/';

% meshname      = 'pensatore';
% ffield_folder = '../ioq_benchmark/results/ashish_nob_02/ffields/';

%meshname      = 'fertility_tri';
%ffield_folder = '../ioq_benchmark/results/ashish_nob_02/ffields/';

N = 4;
%rotate = pi / 4;
rotate = [];

%alg           = 'IOQ_conn_1a';
%alg           = 'JL_IOQ_eps0.5';
%alg            = 'GO';
%alg           = 'MIQ';
algs = {'JL_IOQ_eps0.5'};%, 'GO'};
alg_fname = {'IOQe0_5'};%,'GO'};
out_folder    = 'results/quad/lores/';
if COMPUTE    
    for i = 1:length(algs)
        alg = algs{i};
        name          = [meshname '_' alg];

        grad_size     = 30;
        %qexbin        = '../../ext/libQEx/bin_win64/QEX_bin.exe';
        %QEXBIN = fullfile('.', 'ext', 'libQEx', 'bin_win64', 'QEX_bin.exe');
        qexbin = fullfile('..', '..', 'ext', 'libQEx', 'bin_win64', 'QEX_bin.exe');

        in_ffield = fullfile(ffield_folder, [name '.ffield']);
        disp(['Loading ', in_ffield, '...'])
        m = Mesh(in_ffield); nf = m.nF; nv = m.nV;

        if isempty(rotate)
            R = m.ffield_vectors(1:nf, :);
        else
            theta = m.ffield_angles + rotate;
            R = angles_to_ffield(theta, m.local_frames, N);
            R = R(1:nf, :);
        end

        disp('Running MIQ param...')
        [uv, fuv] = MIQ_param_mex_bin(m.V, m.F, R, grad_size);

        disp('Saving MIQ param...')
        %options.nm_file = ['ext/rosy2gridparam' 'cross.png'];
        options.nm_file = './cross.png';
        options.face_texcorrd = fuv;
        options.object_texture = uv;
        write_obj_for_quad(out_folder, ...
            [name '_gridparam.obj'], ...
            m.V, ...
            m.F, ...
            options);

        disp('Creating quad mesh...')
        param_fp = fullfile(out_folder, [name, '_gridparam.obj']);
        out_quad = fullfile(out_folder, [name, '_quad.obj']);
        system([qexbin, ' ', param_fp, ' ', out_quad]);

        disp('Done.')
    end
end

% % Save camera
% name = [meshname '_MIQ'];
% out_quad = fullfile(out_folder, [name, '_quad']);
% m = MESHQ(out_quad);
% figure; MESH_VIS.mesh( m );
% camf = MESH_VIS.get_camera( gca );
% save([meshname '_camf.mat'],'camf');
% MESH_VIS.set_camera( gca, camf );
% camb = MESH_VIS.get_camera( gca );
% save([meshname '_camb.mat'],'camb');
% camf2 = MESH_VIS.get_camera( gca );
% save([meshname '_camf2.mat'],'camf2');
% MESH_VIS.set_camera( gca, camf2 );
% % 
cams = {};
load( [meshname '_camf.mat'] );
cams{1} = camf;
load( [meshname '_camb.mat'] );
cams{2} = camb;
% load( [meshname '_camf2.mat'] );
% cams{3} = camf2;

colors = [...
    .6,.8,1;...    
    .72,.68,1;...
    1,.9,.7;...
    1,1,1';...
    ];

out_folder    = 'results/quad/lores/';
pannames = {};
pannamesw = {};
if VIS
    for i = 1:length(algs)
        alg = algs{i};
        name = [meshname '_' alg];
        out_quad = fullfile(out_folder, [name, '_quad']);
        m = MESHQ(out_quad);

        outnames = {}; 
        outnamesw = {};
        for j = 1:length(cams)
            outname = fullfile(out_folder, [meshname '_' alg_fname{i}...
                '_quad_cam' num2str(j)]);
            MESH_IO.wfigs(outname,m,colors(i,:),...
                'Plot','mesh','Montage',0,...
                 'OpenGL',1,'Camera',cams{j});        
              
            outnames{end+1} = [outname '.png'];
        end
        
        for j = 1:length(cams)
            outname = fullfile(out_folder, [meshname '_' alg_fname{i}...
                '_quad_white_cam' num2str(j)]);
            MESH_IO.wfigs(outname,m,[1,1,1],...
                'Plot','mesh','Montage',0,...
                 'OpenGL',0,'Camera',cams{j});        
              
            outnamesw{end+1} = [outname '.png'];
        end
        
        out_pan = fullfile(out_folder, [meshname '_' alg_fname{i} '.jpg']);
        pannames{end+1} = out_pan;
        cmd = ['"C:/Program Files (x86)/IrfanView/i_view32.exe" /panorama=(1'];
        for j=1:length(cams)
            cmd = [cmd ',' outnames{j}];
        end
        cmd = [cmd ') /resample /convert=' out_pan];
        system(cmd)        

        out_panw = fullfile(out_folder, [meshname '_white_' alg_fname{i} '.jpg']);
        pannamesw{end+1} = out_panw;
        cmd = ['"C:/Program Files (x86)/IrfanView/i_view32.exe" /panorama=(1'];
        for j=1:length(cams)
            cmd = [cmd ',' outnamesw{j}];
        end
        cmd = [cmd ') /resample /convert=' out_panw];
        system(cmd)        
    end
end

out_pan = fullfile(out_folder, [meshname '.jpg']);
cmd = ['"C:/Program Files (x86)/IrfanView/i_view32.exe" /panorama=(1'];
for i=1:length(algs)
    cmd = [cmd ',' pannames{i}];
end
cmd = [cmd ') /convert=' out_pan];
system(cmd)

out_panw = fullfile(out_folder, [meshname '_w.jpg']);
cmd = ['"C:/Program Files (x86)/IrfanView/i_view32.exe" /panorama=(1'];
for i=1:length(algs)
    cmd = [cmd ',' pannamesw{i}];
end
cmd = [cmd ') /convert=' out_panw];
system(cmd)

