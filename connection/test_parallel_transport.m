m = Mesh('../data/torus_s0.off');
nv = m.nV; nf = m.nF; ne = m.nE;
EF = m.EFAdj;
FE = m.FEAdj;
EV = m.EVAdj;
FF = m.FFAdj;
vf_1ring = m.vf_1ring;

[local_frames, frame_diffs] = create_local_frames(m);
Ad_vert = get_gaussian_curvature(m);

ve_1ring = cell(nv, 1);
for e = 1:ne
    v1 = EV(e, 1);
    v2 = EV(e, 2);
    ve_1ring{v1}(end+1) = e;
    ve_1ring{v2}(end+1) = e;
end

FF_to_frame_diff = sparse(m.nF, m.nF);
be = m.boundaryEdge;
for eid = 1:m.nE
    if be(eid)
        continue;
    end
    fi = EF(eid, 1);
    fj = EF(eid, 2);
    FF_to_frame_diff(fi, fj) = frame_diffs(eid);
    FF_to_frame_diff(fj, fi) = -frame_diffs(eid);
end

for v = 1:nv
    %v = 111;
	e = ve_1ring{v}(1);
	f1 = EF(e, 1);
	f2 = EF(e, 2);
	theta1 = 0;
	theta2 = theta1 + FF_to_frame_diff(f1, f2);
	vec1 = local_angle_to_gvector(f1, theta1, local_frames);
	%vec2 = local_angle_to_gvector(f2, theta2, local_frames);
	visited = [f1];


	while true
		adj_faces = FF(f2, :);
		f_next = setdiff(intersect(adj_faces, vf_1ring{v}), visited);
		if isempty(f_next)
			f_next = f1;
            theta2 = theta2 + FF_to_frame_diff(f2, f_next);
            break;
		end
		theta2 = theta2 + FF_to_frame_diff(f2, f_next);

        visited(end+1) = f2;
		f2 = f_next;
	end

	%res1(v) = mod(theta2 - theta1, 2*pi);
    res1(v) = -theta2;
    vec2 = local_angle_to_gvector(f1, theta2, local_frames);
    res2(v) = vec_vec_angle(vec1, vec2);
    res2(v) = res2(v) * sign(dot(cross(vec2, vec1), m.FNormals(f1, :)));
    %res1(v)
    %res2(v)
    %Ad_vert(v)
end

norm(res1 - Ad_vert)
norm(res2 - Ad_vert)

compare = [Ad_vert, res1(:), res2(:)];
