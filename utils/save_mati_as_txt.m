function save_mati_as_txt(A, filename)
    %function save_mati_as_txt(A, filename)
    % Save the given (integer) matrix as a txt file.

    fid = fopen(filename,'wt');

    for ii = 1:size(A,1)
        fprintf(fid,'%g ', A(ii,:));
        fprintf(fid,'\n');
    end

    fclose(fid);

end

