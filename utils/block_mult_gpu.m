function [result] = block_mult_gpu(A, B, blksize, type)
%BLOCK_MULT_GPU Calculate A*B on the gpu in blocks to reduce memory
%consumption.

if nargin < 4
    type = 'double';
end

[n, m] = size(A);
[r, c] = size(B);
result = zeros(n, c, type);

A = gpuArray(A);

if blksize >= c
    result = gather(A * B);
    return;
end

tmp = gpuArray(zeros(r, blksize, type));
for i = 0:(c/blksize - 1)
    columns = i*blksize + (1:blksize);
    tmp = B(:, columns);
    result(:, columns) = gather(A * tmp);
end

% if l is not a multiple of blksize
columns = (columns(end)+1):c;
if ~isempty(columns)
    tmp(:, columns) = B(:, columns);
    result(:, columns) = gather(A * tmp(:, columns));
end

assert(size(result, 1) == n)
assert(size(result, 2) == c)

end

