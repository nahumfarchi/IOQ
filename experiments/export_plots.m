%%
BASE_FOLDER = 'imgs';
if ~exist(BASE_FOLDER, 'dir')
    mkdir(BASE_FOLDER);
end

filepaths = get_filepaths('.', '.png');

n_files = numel(filepaths);
for i = 1:n_files
    fp = filepaths{i};
    [fp_folder, name, ext] = fileparts(fp);
    if contains(fp_folder, BASE_FOLDER)
        continue
    end
    out_folder = fullfile(BASE_FOLDER, fp_folder);
    if ~exist(out_folder, 'dir')
        mkdir(out_folder);
    end
    fp_out = fullfile(out_folder, [name, ext]);
    if copyfile(fp, fp_out) == 0
        fprintf('Failed to copy %s\n', name);
    end
    
end

%%


