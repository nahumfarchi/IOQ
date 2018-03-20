filepaths = get_filepaths('../data/simple', 'obj');
for i = 1:numel(filepaths)
    ifp = filepaths{i};
    [p, n, e] = fileparts(ifp);
    ofp = fullfile(p, sprintf('%s.off', n));
    if exist(ofp, 'file')
        continue
    end
    
    m = Mesh();
    m.loadTM(ifp);
    m.saveTM(ofp);
end