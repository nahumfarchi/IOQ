function [X,T] = read_off2(filename)

fid = fopen(filename,'r');
if( fid==-1 )
    error('Cannot open the file.');
    return;
end

str = fgets(fid);   % -1 if eof
if ~strcmp(str(1:3), 'OFF')
    error('The file is not a valid OFF one.');    
end

str = fgets(fid);
sizes = sscanf(str, '%d %d', 2);
while length(sizes) ~= 2
    str = fgets(fid);
    sizes = sscanf(str, '%d %d', 2);
end
nv = sizes(1);
nf = sizes(2);

% Read vertices
[X,cnt] = fscanf(fid,'%lf %lf %lf\n', [3,nv]);
if cnt~=3*nv
    warning('Problem in reading vertices.');
end
X = X';

[T,cnt] = fscanf(fid,'3 %ld %ld %ld\n', [3,inf]);
T = T'+1;
T = double(T);

fclose(fid);