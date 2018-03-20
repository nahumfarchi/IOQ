function [latex, T] = bold_latex_table(row_titles, col_titles, X, opt)
%BOLD_LATEX_TABLE Summary of this function goes here
%   Detailed explanation goes here

n = size(X, 1);
m = size(X, 2);
T = cell(n+1, m+1);
T(1, 2:end) = col_titles;
T(2:end, 1) = row_titles;
T(2:end, 2:end) = arrayfun(@(x) num2str(x), X, 'UniformOutput', false);
if isfield(opt, 'min') && opt.min
    [~, bold] = min(X, [], 2);
elseif isfield(opt, 'max') && opt.max
    [~, bold] = max(X, [], 2);
else
    error('Either opt.min or opt.max should be true');
end
for i = 1:n
    j = bold(i);
    T{i+1, j+1} = ['$\textbf{', T{i+1,j+1}, '}$'];
end

opt.data = T;
latex = latexTable(opt);

end

