function [] = draw_local_frames(mesh, local_frames, scale)
    % function [] = draw_local_frames(mesh, local_frames, scale)

    if nargin < 3
        scale = mesh.avg_length;
    end

    mesh.drawFaceField(local_frames(1:mesh.nF, :), 'AutoScale', 'on', 'AutoScaleFactor', scale)
    hold on
    mesh.drawFaceField(local_frames(mesh.nF+1:end, :), 'AutoScale', 'on', 'AutoScaleFactor', scale)
    hold off

end

