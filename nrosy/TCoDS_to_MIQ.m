function [theta, p] = TCoDS_to_MIQ(x, k, r, theta0, Kg, d0, d1, p_fixed)
%function [theta, p] = TCoDS_to_MIQ(x, k, r, theta0, Kg, d0, d1, p_fixed)

nE = size(d0, 1);
nF = size(d1, 1);
pf_inds = ~p_fixed; % free periods
pc_inds = ~pf_inds; % fixed periods
EPS = 1e-10;

% Solve for p
b = (2/pi)*(Kg-d0'*r) - k;
p = zeros(nE, 1);
p(pf_inds) = d0(pf_inds,:)' \ b;

check_norm('d0''*p', 'b');
%disp(norm(d0'*p -b))
%assert(norm(d0'*p -b) < EPS)

% Solve for theta
tf_inds = 2:nF;
tc_inds = 1;
theta = zeros(nF, 1);
theta(tc_inds) = theta0;
theta(tf_inds) = d1(tf_inds,:)' \ (-x-r-(pi/2)*p - d1(tc_inds,:)'*theta(tc_inds));

check_norm('d1''*theta', '(-x-r-(pi/2)*p)');

%disp(['norm(d1^T*res4.thetas - (-y-r-(pi/2)*p)) = ', num2str(norm(d1'*res4.thetas - (-y-r-(pi/2)*p)))])
%disp(['norm(res1.thetas - res4.thetas) = ', num2str(norm(res1.thetas - res4.thetas))])

% Extract p by backslash
% b = (2/pi)*(Kg-d0'*r) - k;
% res4.periods = zeros(m.nE, 1);
% res4.periods(pf_inds) = d0(pf_inds,:)' \ b;
% 
% assert(norm(d0'*res4.periods - b) < EPS)
% disp(['norm(d0^T*res4.periods - b) = ', num2str((norm(d0'*res4.periods - b)))])
% assert(norm(res1.periods-res4.periods) < 1e-6)
% disp(['norm(res1.periods-res4.periods)', num2str(norm(res1.periods-res4.periods))])
% 
% tf_inds = 2:m.nF;
% tc_inds = 1;
% res4.thetas = zeros(m.nF, 1);
% res4.thetas(tc_inds) = theta0;
% res4.thetas(tf_inds) = d1(tf_inds,:)' \ (-y-r-(pi/2)*res4.periods-d1(tc_inds,:)'*res4.thetas(tc_inds));
% 
% assert(norm(d1'*res4.thetas - (-y-r-(pi/2)*p)) < 1e-6)
% assert(norm(res1.thetas - res4.thetas) < 1e-5)
% disp(['norm(d1^T*res4.thetas - (-y-r-(pi/2)*p)) = ', num2str(norm(d1'*res4.thetas - (-y-r-(pi/2)*p)))])
% disp(['norm(res1.thetas - res4.thetas) = ', num2str(norm(res1.thetas - res4.thetas))])
% 
% res4.E = E_MIQ(m, res4.thetas, r, res4.periods, N);
% res4.nrosy = angles_to_nrosy(res4.thetas, res1.local_frames, N);
% 
% disp(['norm(res4.nrosy - res1.nrosy) = ', num2str(norm(res4.nrosy-res1.nrosy))])

end

