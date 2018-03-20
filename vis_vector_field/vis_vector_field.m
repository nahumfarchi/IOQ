function [status, result] = vis_vector_field(mesh_fp, R, S, varargin)
%function [status, result] = vis_vector_field(mesh_fp, rep_fp, ab_fp, varargin)
%
% mesh_fp - path to .off file
% R  - nf x 3 field representative vectors 
% S  - nv x 1 singularity vector (optional)
%
% Examples:
%   vis_vector_field('pensatore.off', R, S);
%   vis_vector_field('pensatore.off', ...
%       R, ...
%       S, ...
%       'degree', 3, ...
%       'lw', 2, ...
%       'max_anim_t', 10, ...
%       'percentage', 0.01);

if nargin < 3
    S = [];
end

p = inputParser;
addOptional(p, 'degree', 4);       % field's degree
addOptional(p, 'lw', 0.5);         % linewidth
addOptional(p, 'miq', false);      % run miq to generate the field
addOptional(p, 'max_anim_t', 20);  % number of streamline iterations
addOptional(p, 'percentage', 0.1); % streamlines sparsity
addOptional(p, 'out', []);         % output png file
addOptional(p, 'width', 5000);     % png width
addOptional(p, 'CM', []);          % colormap matrix
addOptional(p, 'n_colors', 32);    % number of colors to use from the colormap
addOptional(p, 'mesh_color', []);  %
addOptional(p, 'pos_color', []);   % positive singularity color
addOptional(p, 'neg_color', []);   % negative singularity color
addOptional(p, 'sing_size', 30);
addOptional(p, 'angle', []);       % model angle
addOptional(p, 'camera_zoom', 1);
addOptional(p, 'Close', false);
addOptional(p, 'wireframe_only', false);
addOptional(p, 'BinaryPath', 'vis_vector_field_bin.exe');
parse(p, varargin{:});
opt = p.Results;

% Save input matrices as .dmat file
[~, meshname, ~] = fileparts(mesh_fp);
tmp_folder = 'tmp';
mkdir(tmp_folder);

nf = size(R, 1);
R_fp = fullfile(tmp_folder, [meshname '_R.dmat']);
fid = fopen(R_fp, 'w+');
fprintf(fid, '%d %d\n', nf, 3);
fprintf(fid, '%.10g %.10g %.10g \n', R');
fclose(fid);

nv = size(S, 1);
S_fp = fullfile(tmp_folder, [meshname '_S.dmat']);
fid = fopen(S_fp, 'w+');
fprintf(fid, '%d %d\n', nv, 1);
fprintf(fid, '%.10g \n', S');
fclose(fid);

CM = opt.CM;
if ~isempty(CM)
    CM_fp = fullfile(tmp_folder, [meshname '_CM.dmat']);
    fid = fopen(CM_fp, 'w+');
    fprintf(fid, '%d %d\n', size(CM, 1), 3);
    fprintf(fid, '%.10g %.10g %.10g \n', CM');
    fclose(fid);
end

% Create command string

bin_path = opt.BinaryPath;
cmd_str = [bin_path, ' -m ', mesh_fp, ...
                     ' -r ', R_fp, ...
                     ' -d ', num2str(opt.degree), ...
                     ' --lw ', num2str(opt.lw), ...
                     ' --max_anim_t ', num2str(opt.max_anim_t), ...
                     ' --percentage ', num2str(opt.percentage), ...
                     ' --xw ', num2str(opt.width), ...
                     ' --n_colors ', num2str(opt.n_colors), ...
                     ' --sing_size ', num2str(opt.sing_size), ...
                     ' --camera_zoom ', num2str(opt.camera_zoom)];
if opt.miq
    cmd_str = [cmd_str, ' --miq'];
end
if ~isempty(S)
    cmd_str = [cmd_str, ' -s ', S_fp];
end
if ~isempty(opt.out)
    cmd_str = [cmd_str, ' --out ', opt.out];
end
if ~isempty(opt.CM)
    cmd_str = [cmd_str, ' --cm ', CM_fp];
end
if ~isempty(opt.mesh_color)
    cmd_str = [cmd_str, ' --mesh_color ', opt.mesh_color];
end
if ~isempty(opt.pos_color)
    cmd_str = [cmd_str, ' --pos_color ', opt.pos_color];
end
if ~isempty(opt.neg_color)
    cmd_str = [cmd_str, ' --neg_color ', opt.neg_color];
end
if ~isempty(opt.angle)
    cmd_str = [cmd_str, ' --angle ', opt.angle];
end
if opt.Close
    cmd_str = [cmd_str, ' --close'];
end
if opt.wireframe_only
    cmd = [cmd_str, ' --wireframe_only'];
end

% Run exe

disp(['Running ', cmd_str, '...']) 
%[status, result] = system([cmd_str]);
system([cmd_str])

if ~isempty(opt.out)
    RemoveWhiteSpace([], 'file', opt.out);
end

end