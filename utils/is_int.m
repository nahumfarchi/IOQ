function [res] = is_int(x)
%IS_INT Summary of this function goes here
%   Detailed explanation goes here

res = abs(floor(x)-x) < 1e-10;

end

