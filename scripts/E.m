%function [f, g] = E(x, m, d0, A)
function [f, g] = E(x, A, ba, B, bb)
%E Summary of this function goes here
%   Detailed explanation goes here

%f = norm([d0'*d0, zeros(m.nV, m.nV);...
%          zeros(m.nV, m.nV), eye(m.nV)]*x)^2;
%f = norm([d0'*d0, zeros(m.nV, m.nV); ...
%          zeros(m.nV, m.nV), zeros(m.nV, m.nV)]*x)^2;

% just u
%f = norm(d0'*d0*x)^2;

% u and k, no norm on k
%A = [d0'*d0, zeros(m.nV, m.nV);
%     zeros(m.nV, m.nV), zeros(m.nV)];
%f = norm(A*x)^2;

% u and k, with norm on both
%A = [d0'*d0, zeros(m.nV, m.nV);
%     zeros(m.nV, m.nV), eye(m.nV)];
f = norm(A*x-ba)^2;

if nargout > 1
    %g = A*x;
    % just u
    %g = 2*(d0'*d0)*(d0'*d0)*x;
    
    % u and k, no norm on k
    %B = [2*(d0'*d0)*(d0'*d0), zeros(m.nV);
    %     zeros(m.nV), zeros(m.nV)];
    %g = B*x;
    
    % u and k, norm on both
    %B = [2*(d0'*d0)*(d0'*d0), zeros(m.nV);
    %     zeros(m.nV), 2*eye(m.nV)];
    g = B*x - bb;
end

end

