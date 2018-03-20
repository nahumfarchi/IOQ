experiments_setup;

overwrite = true;
filepaths = get_filepaths(out_folder_ffields, ffield_ext);
for i = 1:numel(filepaths)
    fp = filepaths{i};
    [~, meshname, ~] = fileparts(fp);
    if ~overwrite && ...
            exist([out_folder_gridparams, meshname, obj_ext], 'file')
        disp('Skipping (already exists)...')
        continue
    end
    
    disp(['loading ',  fp, '...'])
    m = Mesh();
    m.loadFField(fp);
    [out_folder, name, ~] = fileparts(fp);
    out_fp = fullfile(out_folder_intsing, [name '.txt']);
    disp(['Saving ', fp_out, '...'])
    m.saveSingInteger(out_fp);
    %m.saveFField(fp);
    
    % disp('Running MIQ param...')
    % R = m.ffield_vectors(1:m.nF, :);
    % [uv, fuv] = MIQ_param_mex_bin(m.V, m.F, R, 75);
    
    % disp('Saving...')
    % %options.nm_file = ['ext/rosy2gridparam' 'cross.png'];
    % options.nm_file = './cross.png';
    % options.face_texcorrd = fuv;
    % options.object_texture = uv;
    
    % write_obj_for_quad(out_folder_gridparams, ...
    %     [meshname, obj_ext], ...
    %     m.V, ...
    %     m.F, ...
    %     options);
    
    % disp('Saving quad...')
    % in_obj = [out_folder_gridparams, meshname, obj_ext];
    % out_obj = [out_folder_quads, meshname, obj_ext];
    % system([QEXBIN, ' ', in_obj, ' ', out_obj]);
end