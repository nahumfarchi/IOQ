function [path] = find_data_folder()
%FIND_DATA_FOLDER Return the folder where the data sits.
if exist(fullfile('.', 'data'), 'dir')
    path = fullfile('.', 'data');
elseif exist(fullfile('..', 'data'), 'dir')
    path = fullfile('..', 'data');
elseif exist(fullfile('..', '..', 'data'), 'dir')
    path = fullfile('..', '..', 'data');
else
    error('Data folder does not exist')
end

end

