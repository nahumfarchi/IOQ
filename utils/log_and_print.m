function [count] = log_and_print(fid, format, varargin)
    %LOG_AND_PRINT Summary of this function goes here
    %   Detailed explanation goes here
    
    count = 0;
    if fid > 0
        count = fprintf(fid, format, varargin{:});
    end
    count = count + fprintf(format, varargin{:});
    
end

