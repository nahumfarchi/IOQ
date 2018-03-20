p = find_data_folder();
FNAME = 'sphere_s0.off';
fp = fullfile(p, FNAME);

N = 4;
%f0 = [1, 3];
%f0 = [1, 7];
f0 = [1, 19];
%v0 = normalize_rows([0, 1, 0; 1, 1, 0]);
tc1 = pi+0.1;
tc2 = pi-0.1;
v0 = [cos(tc1), sin(tc1), 0; cos(tc2), sin(tc2), 0];

m = Mesh();
m.loadTM(fp);
scale_factor = m.avg_length / 6;