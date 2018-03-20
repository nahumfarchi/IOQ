p = find_data_folder();
FNAME = 'squares5_row.off';
fp = fullfile(p, FNAME);

N = 4;
fids_list = [1, 19];

m = Mesh();
m.loadTM(fp);
scale_factor = m.avg_length / 6;
[d0, d1] = get_exterior_derivatives(m);

[~, ~, ~, ~, ~, local_frames, r, ~, ~, ~] = nrosy_wrapper(fp, [1], [1,0,0], N);
%thetas_list = [0.25*pi-0.1, 0.75*pi];
thetas_list = [0.5*pi-1e-5, 0.25*pi];
tmp_thetas = nan(m.nF, 1);
tmp_thetas(fids_list) = thetas_list;
v0 = local_angles_to_gvectors(m, tmp_thetas, local_frames, 1);

[C, bc] = create_constraints_mat(m, fids_list, thetas_list, r);

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
gamma = logspace( -2, 2, 20 );
l2norm = zeros(size(gamma));
l1norm = zeros(size(gamma));
fprintf( 1, '   gamma       norm(x,1)    norm(A*x-b)\n' );
fprintf( 1, '---------------------------------------\n' );
for k = 1:length(gamma),
    fprintf( 1, '%8.4e', gamma(k) );
    cvx_begin
        variable x(m.nE);
        minimize( norm(x(nonboundary_edge_ids), 2) + ...
                  gamma(k)*norm(d0(nonboundary_edge_ids, nonboundary_vert_ids)'*x(nonboundary_edge_ids), 1) )
        subject to
            C*x == bc
    cvx_end
    l1norm(k) = norm(d0(nonboundary_edge_ids, nonboundary_vert_ids)'*x(nonboundary_edge_ids), 1);
    l2norm(k) = norm(x(nonboundary_edge_ids), 2);
    fprintf( 1, '   %8.4e   %8.4e\n', l1norm(k), l2norm(k) );
end
figure(5)
plot( l1norm, l2norm );
xlabel( 'norm(x,1)' );
ylabel( 'norm(A*x-b)' );
grid on