%% Settings
FNAME = 'sphere_s0.off';
%FNAME = 'sphere_s1.off';
%FNAME = 'ellipsoid_s4r.off';
%FNAME = 'round_cuber.off';
%FNAME = 'sphere_s0.off';
%FNAME = 'torus_fat_r2.off';
%FNAME = 'bunny.off';
VERBOSE = true;
EPS = 1e-9;

% Find mesh file
p = find_data_folder();
fp = fullfile(p, FNAME);

% Load mesh
m = Mesh();
m.loadTM(fp);
scale = 2 * m.avg_length;

[d0, d1] = get_exterior_derivatives(m);
%d0 = -d0;
%d1 = -d1;
N = 4;          % Direction field degree (i.e., 4 for cross field)
f0 = [1];       % Starting face
%v0 = [0, 0, 1]; % Starting direction
v0 = m.V(m.F(f0, 2), :) - m.V(m.F(f0, 1), :);

%% libigl MIQ (greedy rounding)
print_header('libigl MIQ')
res1 = nrosy_wrapper(fp, f0, v0, N);

%% TCoDS with libigl singularities
I = [res1.S_inds, -res1.S(res1.S_inds)];
%I = [1, 2];
theta0 = [res1.local_frames(1,:); res1.local_frames(m.nF+2,:)] * v0';
theta0 = atan2(theta0(2), theta0(1));
[res2.x, res2.E] = create_trivial_connection(m, I, N, VERBOSE);
[res2.nrosy, res2.theta] = create_direction_field(m, res2.x, f0, theta0, N, res1.local_frames, res1.frame_diffs); %TODO use same angles as in MIQ above

%% MIQ --> TCoDS
%Kg = get_gaussian_curvature(m);
Kg = res1.Kg;
r = res1.frame_diffs;
p = res1.p;
theta = res1.theta;
% 
% %y = (-1)*r + (-1)*2*pi*p/N - d1'*theta;
% %k0 = (Gk-d0'*r)/(2*pi);
% %k = k0-0.25*d0'*p;
% 
% k0 = (-Gk+d0'*r) / (2*pi);
% % TODO check that it is an integer. Does r have enough percision?
% k = k0 + 0.25*d0'*p;
% y = r + 2*pi*p/N + d1'*theta;
% 
% %assert(norm(d0'*res3.x - (-Gk+2*pi*res1.S)))
% disp(['|2*K_G| = ', num2str(norm(2*Gk))])
% 

% k is not integer
% |d0^Ty - (-K_G+(pi/2)k)| = 1.6147e-14
% res3.E == res1.E
%k = (2/pi)*(Gk+d0'*r)+d0'*p;
%y = d1'*theta + r + (pi/2)*p;

% k is integer
% |d0^Ty - (-K_G+(pi/2)k)| = 1.9775
% res3.E == res1.E
%k = (2/pi)*(-Gk+d0'*r)+d0'*p;
%y = d1'*theta + r + (pi/2)*p;

% k is integer
% |d0^Ty - (-K_G+(pi/2)k)| = 134.1325
% res3.E == res1.E
%k = (2/pi)*(Gk-d0'*r)+d0'*p;
%y = d1'*theta + r + (pi/2)*p;

% k is integer
% |d0^Ty - (-K_G+(pi/2)k)| = 2.0304e-14
% res3.E != res1.E

% works best atm
%k = (2/pi)*(res1.Kg-d0'*r)-d0'*p;
%y = -d1'*theta - r - (pi/2)*p;

[y, k] = MIQ_to_TCoDS(theta, p, r, Kg, d0, d1);

disp(['|d0^Ty - (-K_G+(pi/2)k)| = ', num2str(norm(d0'*y - (-Kg+(pi/2)*k)))])

res3.x = y;
res3.nrosy = create_direction_field(m, res3.x, f0, theta0, N, res1.local_frames, res1.frame_diffs);
res3.E = norm(res3.x)^2;


%% TCoDS --> MIQ
pf_inds = ~res1.period_is_fixed;
pc_inds = ~pf_inds;

%LB = -10*ones(m.nE, 1);
%UB = 10*ones(m.nE, 1);

% LP not feasible
%Aeq = d0(inds,:)';
%beq = (2/pi)*(Kg-Aeq*r(inds)) - k;
%res4.p = zeros(m.nE, 1);

% Extract p by solving a LP
%Aeq = d0';
%beq = (2/pi)*(Kg-Aeq*r) - k;
%% Make sure fixed periods are zero
%D = zeros(m.nE, 1);
%D(pc_inds) = 1;
%D = diag(D);
%Aeq = [Aeq; D];
%beq = [beq; zeros(m.nE, 1)];
%res4.p = linprog(ones(m.nE, 1), [], [], Aeq, beq);

% Extract p by backslash
b = (2/pi)*(Kg-d0'*r) - k;
res4.p = zeros(m.nE, 1);
res4.p(pf_inds) = d0(pf_inds,:)' \ b;

assert(norm(d0'*res4.p - b) < EPS)
disp(['norm(d0^T*res4.p - b) = ', num2str((norm(d0'*res4.p - b)))])
assert(norm(res1.p-res4.p) < 1e-6)
disp(['norm(res1.p-res4.p)', num2str(norm(res1.p-res4.p))])

tf_inds = 2:m.nF;
tc_inds = 1;
res4.theta = zeros(m.nF, 1);
res4.theta(tc_inds) = theta0;
res4.theta(tf_inds) = d1(tf_inds,:)' \ (-y-r-(pi/2)*res4.p-d1(tc_inds,:)'*res4.theta(tc_inds));

assert(norm(d1'*res4.theta - (-y-r-(pi/2)*p)) < 1e-6)
assert(norm(res1.theta - res4.theta) < 1e-5)
disp(['norm(d1^T*res4.theta - (-y-r-(pi/2)*p)) = ', num2str(norm(d1'*res4.theta - (-y-r-(pi/2)*p)))])
disp(['norm(res1.theta - res4.theta) = ', num2str(norm(res1.theta - res4.theta))])

res4.E = E_MIQ(m, res4.theta, r, res4.p, N);
res4.nrosy = angles_to_nrosy(res4.theta, res1.local_frames, N);

disp(['norm(res4.nrosy - res1.nrosy) = ', num2str(norm(res4.nrosy-res1.nrosy))])

%res4.p(inds) = linprog(ones(sum(inds),1), [], [], Aeq, beq);
%res4.p(inds) = d0(inds,:)' \ ((2/pi)*(Kg-d0(inds,:)'*r(inds))); %(4*(k-(Kg-d0(inds,:)'*r(inds))/(2*pi)));
%norm(k - ((2/pi)*(Kg-d0'*r) - d0'*p))
%norm(d0(inds,:)'*res4.p(inds) - (4*(k-(Kg+d0(inds,:)'*r(inds))/(2*pi))))


%res4.p = p;
%inds = 2:m.nF;
%res4.theta = zeros(m.nF, 1);
%res4.theta(1) = theta0;
%res4.theta(inds) = d1(inds,:)' \ (y - r - (pi/2)*eye(length(res4.p))*res4.p - d1(1,:)'*theta(1));
%norm(d1(inds,:)'*res4.theta(inds) - (y-r-(pi/2)*eye(length(res4.p))*res4.p) - d1(1,:)'*theta(1))
%norm(d1'*res4.theta - (y-r-(pi/2)*eye(length(res4.p))*res4.p))

%res4.E = E_MIQ(m, res4.theta, r, res4.p, N);
%res4.nrosy = angles_to_nrosy(res4.theta, res1.local_frames, N);

%% Plots
nPlots = 2;
mPlots = 2;
spi = 1;
figi = 1;
ha = zeros(nPlots*mPlots,1);
scale = 2*m.avg_length;

figure(figi); clf(figi)
ha(1) = subplot(nPlots, mPlots, spi);
spi = spi + 1;
hold on
draw_constraints(m, f0, v0)
draw_direction_field(m, res1.nrosy, res1.S_inds, scale)
hold off
title({['libigl MIQ (greedy rounding)'], ['E = ', num2str(res1.E)]})

ha(2) = subplot(nPlots, mPlots, spi);
spi = spi + 1;
hold on
draw_constraints(m, f0, v0);
draw_direction_field(m, res2.nrosy, res1.S_inds, scale)
hold off
title({'TCoDS with MIQ singularities', ['E = ', num2str(res2.E)]})

ha(3) = subplot(nPlots, mPlots, spi);
spi = spi + 1;
hold on
draw_constraints(m, f0, v0);
draw_direction_field(m, res3.nrosy, res1.S_inds, scale)
hold off
title({'MIQ(p,theta) --> TCoDS(x,k)', ['E = ', num2str(res3.E)]})

ha(4) = subplot(nPlots, mPlots, spi);
spi = spi + 1;
hold on;
draw_constraints(m, f0, v0);
draw_direction_field(m, res4.nrosy, res1.S_inds, scale)
hold off
title({'TCoDS(x,k) --> MIQ(p,theta)', ['E = ', num2str(res4.E)]})

hlink = linkprop(ha, {'CameraPosition', 'CameraUpVector'});
figi = figi + 1;

%% debug orientation
% figure(figi); clf(figi)
% clear ha
% %ha(1) = subplot(221);
% %m.draw()
% %m.drawLabels(true, true, true)
% 
% ha(1) = subplot(121);
% 
% VERT = 52;
% m.labelEdges(r, '%.2f', [0, 0, 0], 'Color', 'blue');
% m.labelVertices(Kg, '%.2f', [0, 0, 0], 'Color', 'red');
% m.labelEdges(full(d0(:,VERT)), '%.2f', [0,0,0.1*m.avg_length], 'Color', 'g')
% hold on
% m.drawFaceField(res1.local_frames(1:m.nF, :), 'AutoScale', 'on', 'AutoScaleFactor', scale/2)
% m.labelEdges(y, '%.4f', [0,0,0.2*m.avg_length], 'Color', 'c')
% 
% m.labelFaces(res1.theta, '%.2f', [0, 0, 0], 'Color', 'm')
% m.labelEdges(p, '%d', [0,0,-0.1*m.avg_length], 'Color', 'k')
% m.draw();
% 
% P = zeros(m.nE, 3);
% U = zeros(m.nE, 3);
% for eid = 1:m.nE
%     vert1 = m.V(m.EVAdj(eid, 1), :);
%     vert2 = m.V(m.EVAdj(eid, 2), :);
%     u = vert2 - vert1;
%     P(eid, :) = (vert1+vert2)/2;
%     U(eid, :) = u;
% end
% quiver3(P(:,1), P(:,2), P(:,3), U(:,1), U(:,2), U(:,3))
% 
% hold off
% axis off
% title({'blue - r', 'red - Gaussian curvature', 'green - d0(:,52)'})
% 
% EDGES = [271, 269];
% ha(2) = subplot(122);
% m.draw();
% m.drawLabels(true, false, true)
% hold on
% m.labelEdges(full(d0(:,VERT)), '%.2f', [0,0,0.1*m.avg_length], 'Color', 'g')
% m.labelFaces(full(-d1(:,271)), '%d', [0,0,0], 'color', 'm')
% m.labelFaces(full(-d1(:,269)), '%d', [0,0,0.1*m.avg_length], 'Color', 'k')
% m.labelFaces(full(-d1(:,268)), '%d', [0,0,0.2*m.avg_length], 'Color', 'g')
% %for k = length(EDGES)
% %    eid = EDGES(k);
% %    m.labelFaces(full(-d1(:,eid)), '%d', [0,0,0.1*(k-1)*m.avg_length], 'color', 'm')
% %end
% hold on
% 
% axis off
% 
% hlink = linkprop(ha, {'CameraPosition', 'CameraUpVector'});

