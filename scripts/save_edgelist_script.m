folder = 'edgelist';
mkdir(folder);

p = find_data_folder();

m = Mesh();
m.loadTM(fullfile(p, 'bumpy.off'));

edgelist = create_dual_edgelist(m);

path = fullfile(folder, 'bumpy_edgelist.txt');
fid = fopen(path, 'w');
fprintf(fid, '%d\n', size(edgelist, 1));
fprintf(fid, '%d %d\n', edgelist');
fclose(fid);