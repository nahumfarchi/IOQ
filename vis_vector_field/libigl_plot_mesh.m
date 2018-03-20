function [status, result] = libigl_plot_mesh(mesh_fp, varargin)


p = inputParser;
addOptional(p, 'lw', 0.5);         % linewidth
addOptional(p, 'out', []);         % output png file
addOptional(p, 'width', 5000);     % png width
addOptional(p, 'mesh_color', '1,1,1');
addOptional(p, 'sing_size', 30);
addOptional(p, 'angle', []);       % model angle
addOptional(p, 'camera_zoom', 1);
addOptional(p, 'Close', false);
addOptional(p, 'BinaryPath', 'vis_vector_field_bin.exe');
parse(p, varargin{:});
opt = p.Results;

% Create command string

bin_path = opt.BinaryPath;
cmd_str = [bin_path, ' -m ', mesh_fp, ...
                     ' --lw ', num2str(opt.lw), ...
                     ' --xw ', num2str(opt.width), ...
                     ' --sing_size ', num2str(opt.sing_size), ...
                     ' --camera_zoom ', num2str(opt.camera_zoom), ...
                     ' --wireframe_only '];
if ~isempty(opt.out)
    cmd_str = [cmd_str, ' --out ', opt.out];
end
if ~isempty(opt.mesh_color)
    cmd_str = [cmd_str, ' --mesh_color ', opt.mesh_color];
end
if ~isempty(opt.angle)
    cmd_str = [cmd_str, ' --angle ', opt.angle];
end
if opt.Close
    cmd_str = [cmd_str, ' --close'];
end

% Run exe

disp(['Running ', cmd_str, '...']) 
%[status, result] = system([cmd_str]);
system([cmd_str])

if ~isempty(opt.out)
    RemoveWhiteSpace([], 'file', opt.out);
end

end