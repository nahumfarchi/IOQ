p = find_data_folder();
%FNAME = 'ellipsoid_s4r.off';
%FNAME = 'round_cuber.off';
FNAME = 'sphere_s2.off';
fp = fullfile(p, FNAME);

N = 4;
fids_list = [1, 19];

m = Mesh();
m.loadTM(fp);
%scale_factor = m.avg_length / 6;
[d0, d1] = get_exterior_derivatives(m);

%[~, ~, ~, ~, ~, local_frames, r, ~, ~, ~] = nrosy_wrapper(fp, [1], [1,0,0], N);
%thetas_list = [0.25*pi-0.1, 0.75*pi];
%thetas_list = [0.5*pi-1e-5, 0.25*pi];
%tmp_thetas = nan(m.nF, 1);
%tmp_thetas(fids_list) = thetas_list;
%v0 = local_angles_to_gvectors(m, tmp_thetas, local_frames, 1);
Gk = get_gaussian_curvature(m);

%[C, bc] = create_constraints_mat(m, fids_list, thetas_list, r);

nonboundary_edge_ids = [];
nonboundary_vert_ids = zeros(m.nV, 1);
for eid = 1:m.nE
    if ~m.isBoundaryEdge(eid)
        nonboundary_edge_ids(end+1) = eid;
    else
        v1 = m.EVAdj(eid, 1);
        v2 = m.EVAdj(eid, 2);
        nonboundary_vert_ids(v1) = 1;
        nonboundary_vert_ids(v2) = 1;
    end
end
nonboundary_vert_ids = find(nonboundary_vert_ids == 0);

%% Paramater search
lambda = logspace( -2, 2, 9 );
%lambda = [0, 1, 100];
l2norm = zeros(size(lambda));
l1norm = zeros(size(lambda));
solutions = zeros(m.nE, length(lambda));
fprintf( 1, '   lambda       norm(x,1)    norm(A*x-b)\n' );
fprintf( 1, '---------------------------------------\n' );
figure(1)
for i = 1:length(lambda),
    fprintf( 1, '%8.4e', lambda(i) );
    cvx_begin
        variable x(m.nE)
        minimize( norm(x, 2) + lambda(i)*norm(d0'*x+Gk, 1) )
        %subject to
        %    C*x == bc
    cvx_end
    l1norm(i) = norm(d0'*x+Gk, 1);
    l2norm(i) = norm(x, 2);
    solutions(:, i) = x;
    fprintf( 1, '   %8.4e   %8.4e\n', l1norm(i), l2norm(i) );
    ha(i) = subplot(3,3,i); m.draw(d0'*x, 'EdgeAlpha', 0.0); title(['\lambda=', num2str(lambda(i))]);
    colorbar
end
hlink = linkprop(ha, {'CameraPosition', 'CameraUpVector'});

%%
figure(1)
for i = 1:length(lambda)
    ha(i) = subplot(3,3,i);
    x = solutions(:,i);
    m.draw(d0'*x, 'EdgeAlpha', 0.0)
    title(['\lambda=', num2str(lambda(i))])
    colorbar
end
hlink = linkprop(ha, {'CameraPosition', 'CameraUpVector'});

figure(2)
plot( l1norm, l2norm, '-rx' );
xlabel( 'norm(d0^Tx+Gk,1)' );
ylabel( 'norm(x, 2)' );
title([FNAME, ' (logscale)'])
grid on

figure(3)
loglog( lambda, l1norm, '-rx', lambda, l2norm, '-bx')
xlabel( 'lambda' )
legend('norm(d0^Tx+Gk,1))', 'norm(x, 2)')
title([FNAME, ' (logscale)'])

figure(4)
for i = 1:length(lambda)
    subplot(3,3,i);
    x = solutions(:,i);
    hist(d0'*x,100)
    %axis([-0.5 0.5 0 3500])
    title(['\lambda=', num2str(lambda(i))]);
end

figure(5)
m.draw(Gk, 'EdgeAlpha', 0.0);
title('Gaussian curvature')