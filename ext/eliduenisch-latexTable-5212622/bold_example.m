n = 4;
A = magic(n);
COL_TITLES = {'col1', 'col2', 'col3', 'col4'};
ROW_TITLES = {'row1', 'row2', 'row3', 'row4'};
T = cell(n+1, n+1);
T(1, 2:end) = COL_TITLES;
T(2:end, 1) = ROW_TITLES;
T(2:end, 2:end) = arrayfun(@(x) num2str(x), A, 'UniformOutput', false);
[~, bold] = min(A, [], 2);
for i = 1:n
    j = bold(i);
    T{i+1, j+1} = ['$\textbf{', T{i+1,j+1}, '}$'];
end

clear input;
input.data = table(T);
input.tableCaption = 'Caption';
input.tablePlacement = 'H';
input.tableColumnAlignment = 'c';
latex = latexTable(input);