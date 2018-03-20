function str = create_header(s, margin_width)
    if nargin < 2
        margin_width = 3;
    end
    total_width = 2*margin_width + 2 + length(s);
    row = repmat('*', 1, total_width);
    margin = repmat('*', 1, margin_width);
    str = [row, '\r\n', sprintf('%s %s %s\r\n', margin, s, margin), row, '\r\n'];
end

