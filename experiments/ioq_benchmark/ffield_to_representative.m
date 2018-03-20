overwrite = false;
in_folder = fullfile('results', 'ashish_nob', 'ffields');
out_folder = fullfile('results', 'ashish_nob', 'representative');
mkdir(out_folder);

paths = get_filepaths(in_folder, '.ffield');
n = numel(paths);
progressbar
for i = 1:n
    fp_in = paths{i};
    disp(fp_in)
    [~, name, ~] = fileparts(fp_in);
    fp_out = fullfile(out_folder, [name '.dmat']);
    
    if exist(fp_out, 'file') && ~overwrite
        disp('skipping...')
        progressbar(i / n)
        continue
    end
    
    fid = fopen(fp_out, 'w+');
    if fid < 0
        error('Could not open file')
    end
    
    m = Mesh(fp_in);
    nf = m.nF;
    R = m.ffield_vectors(1:nf, :);
    
    fprintf(fid, '%d %d\n', nf, 3);
    fprintf(fid, '%.10g %.10g %.10g \n', R');
    
    progressbar(i / n)
end

%%
in_folder = fullfile('results', 'ashish_nob', 'ffields');
out_folder = fullfile('results', 'ashish_nob', 'alphabeta');
paths = get_filepaths(in_folder, '.ffield');
n = numel(paths);
mkdir(out_folder);
for i = 1:n
    fp_in = paths{i};
    disp(fp_in)
    [~, name, ~] = fileparts(fp_in);
    fp_out = fullfile(out_folder, [name '.dmat']);
    fid = fopen(fp_out, 'w+');
    if fid < 0
        error('Could not open file')
    end
    
    m = Mesh(fp_in);
    nf = m.nF; ng = m.genus; ng2 = 2*ng; nv = m.nV;
    alpha = zeros(nv, 1);
    beta = zeros(ng2, 1);
    alpha(m.vert_sing(:, 1)) = m.vert_sing(:, 2);
    beta(m.gen_sing(:, 1)) = m.gen_sing(:, 2);
    
    if any(alpha-round(alpha) > 1e-6)
        alphabeta = [m.degree*alpha; m.degree*beta];
    else
        alphabeta = [alpha; beta];
    end
    
    fprintf(fid, '%d %d\n', nv, 1);
    fprintf(fid, '%.10g \n', alphabeta');
end