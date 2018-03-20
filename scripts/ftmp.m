function [f, g] = ftmp(x, d0, Ad)
%FTMP Summary of this function goes here
%   Detailed explanation goes here
    %A = d0;
    %b = [Ad;
    n = length(x)/2;
    d0d0 = d0'*d0;
    A = [d0d0, (-pi/2)*eye(n)];
    b = Ad;
    u = x(1:n);
    k = x(n+1:end);
    
    f = norm(d0*u)^2 + norm(k)^2 + norm(A*x + b)^2;
    g = 2*A'*(A*x+b);
    g(1:n) = g(1:n) + 2*d0d0*u;
    g(n+1:end) = g(n+1:end) + 2*k;
    
    %f = norm(x)^2;
    %g = 2*x;
    
    %n = length(x)/2;
    %u = x(1:n);
    %k = x(n+1:end);
    %d0d0 = d0'*d0;
    
    %f = norm(d0d0*u+Ad-(pi/2)*k)^2;
    %gu = ((d0d0*u+Ad-(pi/2)*k)'*d0d0)';
    %gk = (-pi/2)*(d0d0*u+Ad-(pi/2));
    %g = [gu; gk];

end

