function [fpaths] = get_filepaths(directory, ext, rec)
    % function [fpaths] = get_filepaths(directory, ext)
    % Return a list of all file paths in the given directory and
    % subdirectories (recursive).
    if nargin < 2
        ext = '';
    elseif ext(1) ~= '*'
        if ext(1) ~= '.'
            ext = strcat('*.', ext);
        else
            ext = strcat('*', ext);
        end
    end
    if nargin < 3
        rec = true;
    end
    
    if rec
        [p0, n, e] = fileparts(directory);
        p = strread(genpath(directory), '%s', 'delimiter', ';');
        fpaths = {};
        idx = 1;
        for i = 1:length(p)
            f = dir(fullfile(p{i}, ext));
            for j = 1:length(f)
                fp = fullfile(p{i}, f(j).name);           
                if ~isdir(fp)
                    fpaths{idx} = fp;
                    idx = idx + 1;
                end
            end
        end
    else
        f = dir(fullfile(directory, ext));
        for j = 1:length(f)
            fp = fullfile(directory, f(j).name);
            if ~isdir(fp)
                fpaths{j} = fp;
            end
        end
    end
end

