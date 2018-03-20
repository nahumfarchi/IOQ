%% Setup
global LOG;
LOG = -1;

if ~exist('experiment_name', 'var')
    disp('Set variable experiment_name!')
    return
end

print_header(['experiment_name : ', experiment_name])

%experiment_name = 'bunnies_connectivity_lap';
%data_folder = '../data/genus0';
base_folder = fullfile('..', 'results', 'experiments', experiment_name);
out_folder_ffields = fullfile(base_folder, 'ffields');
face_based = true;
if face_based
    out_folder_nrosy = fullfile(base_folder, 'nrosy_face');
else
    out_folder_nrosy = fullfile(base_folder, 'nrosy_vertex');
end
out_folder_plots = fullfile(base_folder, 'plots');
ext_in = '.ffield';
ext_out = '.rosy';

disp(['Results are being saved into ', base_folder])

if ~exist(out_folder_ffields, 'dir')
    mkdir(out_folder_ffields);
end
if ~exist(out_folder_nrosy, 'dir')
    mkdir(out_folder_nrosy);
end
if ~exist(out_folder_plots, 'dir')
    mkdir(out_folder_plots);
end

LOG = fopen(fullfile(out_folder_ffields, 'log.txt'), 'w');

ffield_filepaths = get_filepaths(out_folder_ffields, ext_in);
n_files = numel(ffield_filepaths);
for r = 1:n_files
    fp_in = ffield_filepaths{r};
    [~, name, ext_in] = fileparts(fp_in);
    fp_out = fullfile(out_folder_nrosy, [name, ext_out]);
    m = Mesh();
    m.loadFField(fp_in);
    m.saveNRosy(fp_out, face_based);
end