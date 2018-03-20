ioq_benchmark_setup;
%experiments_setup;

filepaths = get_filepaths(out_folder_gridparams, obj_ext);
for i = 1:numel(filepaths)
    fp = filepaths{i};
    disp(['loading',  fp, '...'])
    
    %disp('Running MIQ param...')
    %R = m.ffield_vectors(1:m.nF, :);
    %[uv, fuv] = MIQ_param_mex_bin(m.V, m.F, R, 75);
    disp('Converting to quad...')
    [~, meshname, ~] = fileparts(fp);
    out_fp = [out_folder_quads, meshname, obj_ext];
    system([qexbin, ' ', fp, ' ', out_fp]);
    
    %disp('Saving...')
    %options.nm_file = ['ext/rosy2gridparam' 'cross.png'];
    %options.nm_file = './cross.png';
    %options.face_texcorrd = fuv;
    %options.object_texture = uv;
    
    %write_obj(out_folder_quads, [meshname, '.obj'], m.V, m.F, options);
end