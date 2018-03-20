function [x, k] = MIQ_to_TCoDS(theta, p, r, Kg, d0, d1, degree)
% function [x, k] = MIQ_to_TCoDS(theta, p, r, Kg, d0, d1, degree)

if nargin < 7
    degree = 4;
end

% Works with libigl frames
%k = (2/pi)*(Kg-d0'*r) - d0'*p;
%%k(abs(k) < 1e-5) = 0;
%x = -d1'*theta - r - (pi/2)*p;

k = (degree/(2*pi))*(Kg-d0'*r) - d0'*p;
%k(abs(k) < 1e-5) = 0;
x = d1'*theta + r + (2*pi/degree)*p;
x = -x;

end

