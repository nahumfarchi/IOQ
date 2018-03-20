function [f, g] = E_2017_07_16_1105(x, d0, Ad, la)
%E_2017_07_16_1105 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 4
    la = 1;
end

n = length(x) / 2;
u = x(1:n);
k = x(n+1:end);

f = norm(d0*u)^2 + norm(k)^2 + la*norm((d0'*d0)*u+Ad-(pi/2)*k)^2;

d0d0 = d0'*d0;
gu = 2*d0d0*u + ...
        2*la*(u'*d0d0*d0d0*u+d0d0*Ad-(pi/2)*d0d0*k);
gk = 2*k + ...
     pi*la*((pi/2)*k-d0d0*u-Ad);
g = [gu; gk];

%gradu = 2 * ( (d0d0 - pi*la*d0d0)*u + ...
%          (eye(n) - (pi*la)/2 * d0d0 + (pi^2)*la/2 * eye(n))*k + ...
%          la*u'*(d0d0*d0d0)*u + ...
%          la*d0d0*Ad - pi*la*Ad );

end

