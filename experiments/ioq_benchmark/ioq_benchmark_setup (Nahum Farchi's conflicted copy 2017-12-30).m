%% ========================================================================
%  Setup 
%  ========================================================================

global LOG;
global EPS;
global VERBOSE;
global PLOT;
global SAVE;
global DEGREE;
global ASPECT_RATIO;
global RESOLUTION;
VERBOSE = true;
EPS = 1e-9;
PLOT = true;
SAVE = true;
DEGREE = 4;
ASPECT_RATIO = 1;
RESOLUTION = 1024;    
LOG = -1;

F0 = [1];
THETA0 = [0];
% Constraint vector in global coordinates
G_CONSTRAINT_VEC = [1, 0, 0];
N_ITER = 1000;
%LB = 8000; % with 4gb gpu
LB = 20000; % with 12gb gpu

ALG_OPTS = {'ioq_gpu_inv_iter_gpu', ...
        'ioq_block_inv_iter_gpu', ...
        'ioq_cpu_inv_iter_gpu', ...
        'ioq_block_inv_iter_cpu', ...
        'ioq_cpu_inv_iter_cpu', ...
        'miq', ...
        'jl_ioq_iter_gpu', ...
        'jl_ioq_iter_cpu'};
%ALG_OPTS = {'ioq_genus0', ...
%    'ioq_opt1a', ...
%    'ioq_optz2', ...
%    'miq'
ALGS = {'IOQ_conn_opt1a', 'IOQ_conn_opt2z', ...
        'IOQ_cot_opt1a', 'IOQ_cot_opt2z', ...
        'MIQ'};
SYMBOLS = containers.Map(ALG_OPTS, {'s', 'o', 'v', '^', '>', 'd', '+', '*'});
%COLORS = containers.Map(ALGS_OPT, {'r', 'g', 'b', 'm', 'c', 'k'});
ALG_TO_COLOR = containers.Map(ALGS, {'r', 'b', 'k', 'm', 'g'});


%experiment_name = 'bunnies_connectivity_lap';
%data_folder = '../data/bunnies';

%experiment_name = 'bunnies_cot_lap';
%data_folder = '../data/bunnies';

%experiment_name = 'genus0_connectivity_lap';
%data_folder = '../data/genus0';

%experiment_name = 'bunnies';
%data_folder = '../data/bunnies';

experiment_name = 'ashish_nob_small';
%experiment_name = 'ashish_nob';
%experiment_name = 'bunnies';
%experiment_name = 'ashish_nob_small';
%data_folder = '../../../data/ashish_nob/';
data_folder = fullfile('..', '..', '..', 'data', experiment_name);

data_ext = '.off';

%base_folder = fullfile('..', 'results', 'experiments', experiment_name);
base_folder = fullfile('.', 'results', experiment_name);
out_folder_ffields = fullfile(base_folder, 'ffields/');
out_folder_nrosy = fullfile(base_folder, 'nrosy/');
out_folder_plots = fullfile(base_folder, 'plots/');
out_folder_gridparams = fullfile(base_folder, 'gridparams/');
out_folder_quads = fullfile(base_folder, 'quads/');
out_folder_intsing = fullfile(base_folder, 'intsing/');
out_folder_seamless = fullfile(base_folder, 'seamless/');
global OUT_FOLDER_LP;
OUT_FOLDER_LP = fullfile('D:\nahum\laplacian_pinvs\');
folders = {base_folder, ...
           out_folder_ffields, ...
           out_folder_nrosy, ...
           out_folder_plots, ...
           out_folder_gridparams, ...
           out_folder_quads, ...
           out_folder_intsing, ...
           out_folder_seamless, ...
           OUT_FOLDER_LP};

ffield_ext = '.ffield';
nrosy_ext = '.rosy';
obj_ext = '.obj';

% Create folders
for i = 1:numel(folders)
    fld = folders{i};
    if ~exist(fld, 'dir')
        try
            mkdir(fld);
        catch
            disp(['Could not create folder ', fld])
        end
    end
end

LOG = fopen(fullfile(base_folder, 'log.txt'), 'w');

QEXBIN = fullfile('.', 'ext', 'libQEx', 'bin_win64', 'QEX_bin.exe');
% Usage:
%system([QEXBIN, ' ', in_obj_path, ' ', out_obj_path]);

MIQ_COL = 2;