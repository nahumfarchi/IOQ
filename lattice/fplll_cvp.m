function [status, result] = fplll_cvp(M, t)
%FPLLL_CVP Use fplll.exe to find the closest vector point to t in the
%lattice defined by M.
    global FPLLL;
    global TMP_FOLDER;

    if ~exist(TMP_FOLDER, 'dir')
        mkdir(TMP_FOLDER);
    end
    
    fp = fullfile(TMP_FOLDER, 'fplll_cvp_tmp_matrix.txt');
    [fid, msg] = fopen(fp, 'w');
    if fid < 0
        disp(['Failed to open file ', fp])
        disp(msg)
        y = [];
        return;
    end
    
    %fprintf(fid, '[');
    %fprintf(fid, '[%g ]\n', M(:,1:end-1));
    %fprintf(fid, '[%g %g %g]]\n[', M(:, end));
    %fprintf(fid, '%g ', t(1:end-1));
    %fprintf(fid, '%g]', t(end));
    %fprintf(fid, ']\n');
%     M = full(M');
%     fprintf(fid, '[[');
%     fprintf(fid, [repmat('%g ', 1, size(M, 2)) ']\n['], M(:, 1:end-1));
%     fprintf(fid, '%g ', M(:, end));
%     %fprintf(fid, ']]\n[');
%     fprintf(fid, ']]\n[');
%     fprintf(fid, '%g ', t(:));
%     fprintf(fid, ']\n');
    M = round(10000*full(M));
    t = round(10000*t);
    fprintf(fid, '[');
    for i = 1:size(M, 2)
        fprintf(fid, '[');
        fprintf(fid, '%g ', M(:, i));
        if i==size(M, 2)
            fprintf(fid, ']]\n');
        else
            fprintf(fid, ']\n');
        end
    end
    fprintf(fid, '[');
    fprintf(fid, '%g ', t(:));
    fprintf(fid, ']\n');
    
    cmd = [FPLLL, ' ', fp, ' -a cvp'];
    [status, result] = system(cmd);
end

