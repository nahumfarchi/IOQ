function g_vec = local_angle_to_gvector(fid, theta, local_frames)
	nf = size(local_frames, 1) / 2;
	frame = [local_frames(fid, :); ...
			 local_frames(fid+nf, :)];
	g_vec = [cos(theta), sin(theta)] * frame;
end