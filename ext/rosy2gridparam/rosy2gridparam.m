
clearvars; close all;

datadir = './';
outdir = './';

meshname = 'sphere_s2';
[X,T] = read_off([datadir meshname '.off']);
nF = size(T, 1);

% output is saved in outdir with names: meshname '_uv' or '_fuv'
%miqbin = 'MIQ_bin';
%eval(['!' miqbin ' ' datadir ' ' meshname ' ' outdir]);

% load UV and FUV
uv = load([outdir meshname '_uv.txt']);         % size 2d_nv x 2
fuv = load([outdir meshname '_fuv.txt']) + 1;   % size nf x 3

shareddir = './';
options.nm_file = [shareddir 'cross.png'];
options.face_texcorrd = fuv;
options.object_texture = uv;
write_obj( outdir, [meshname '_exe.obj'],X, T, options );


fid = fopen([datadir meshname '_nrosy.txt']);
[R, Rcount] = fscanf(fid, '%g %g %g', [3 nF]); 
R = R';

[uv2, fuv2] = MIQ_param_mex_bin(X, T, R, 75);

check_norm('uv', 'uv2')
check_norm('fuv', 'fuv2')

options.nm_file = [shareddir 'cross.png'];
options.face_texcorrd = fuv2;
options.object_texture = uv2;
write_obj( outdir, [meshname '_mex.obj'],X, T, options );

figure
subplot(121)
patch('Faces',fuv,'Vertices',uv,'FaceColor','w','edgecolor','k'); axis equal
title('MIQ exe')
subplot(122)
patch('Faces',fuv2,'Vertices',uv2,'FaceColor','w','edgecolor','k'); axis equal
title('MIQ mex')