function [folder_path] = create_time_stamped_folder(root, folder, date, hour)
    %function [folder] = create_time_stamped_folder(root)
    
    if nargin < 2
        folder = '';
    end
    if nargin < 3
        date = true;
    end
    if nargin < 4
        hour = true;
    end

    T = datetime;
    if date
        folder_path = fullfile(root, ...
            [num2str(T.Year), '.', sprintf('%02d', month(T)), '.', sprintf('%02d', T.Day)]);
    end
    if hour
        folder_path = fullfile(folder_path, ...
            [sprintf('%02d', T.Hour), sprintf('%02d', (T.Minute))]);
    end
    folder_path = fullfile(folder_path, folder);
    %folder = fullfile(root, ...
    %                      [num2str(T.Year), '.', sprintf('%02d', month(T)), '.', sprintf('%02d', T.Day)], ...
    %                      [sprintf('%02d', T.Hour), sprintf('%02d', (T.Minute))]);
    if ~exist(folder_path, 'dir')
        mkdir(folder_path);
    end

end

