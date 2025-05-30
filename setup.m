%restoredefaultpath
%savepath

global LIBIGL_NROSY; 
global LIBIGL_COMISO;
global TMP_FOLDER; 
global FIG;
global FPLLL;
LIBIGL_NROSY = fullfile('..', 'bin', 'NRoSy_wrapper_bin.exe');
LIBIGL_COMISO = fullfile('libigl', 'build', 'CoMISo_wrapper_bin.exe');
FPLLL = fullfile('..', 'bin', 'fplll', 'fplll.exe');
TMP_FOLDER = fullfile('tmp');

FOLDERS = {'mesh', ...
    'nrosy', ...
    'utils', ...
    'scripts', ...
    'connection', ...
    'libigl', ...
    'phase_unwrapping', ...
    'lattice', ...
    'GODF', ...
    'ext', ...
    'tests', ...
    'IOQ', ...
    'vis_vector_field'};

for i= 1:length(FOLDERS)
    f = FOLDERS{i};
    p = genpath(fullfile(pwd, f));
    fprintf('Adding ''%s'' to path...\n', p);
    addpath(p);
end

rmpath(genpath(fullfile('.', 'ext', 'matlab-toolboxes-master')));

fprintf('Done.\n\n');