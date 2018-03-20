function print_header(s, margin_width)
    if nargin < 2
        margin_width = 3;
    end
    total_width = 2*margin_width + 2 + length(s);
    row = repmat('*', 1, total_width);
    margin = repmat('*', 1, margin_width);
    disp(row);
    fprintf('%s %s %s\n', margin, s, margin);
    disp(row);    
end

