function cam = get_camera(ca, verbose)
    if nargin < 2
        verbose = false;
    end
    cam.pba = get(ca, 'PlotBoxAspectRatio');
    cam.dar = get(ca, 'DataAspectRatio');
    cam.cva = get(ca, 'CameraViewAngle');
    cam.cuv = get(ca, 'CameraUpVector');
    cam.ct = get(ca, 'CameraTarget');
    cam.cp = get(ca, 'CameraPosition');
    
    if verbose
        %disp(['cam.pba = ', cam.pba])
        fprintf('cam.pba = ['); fprintf('%f ', cam.pba);
        fprintf('];\ncam.dar = ['); fprintf('%d ', cam.dar);
        fprintf('];\ncam.cva = ['); fprintf('%f ', cam.cva);
        fprintf('];\ncam.cuv = ['); fprintf('%f ', cam.cuv);
        fprintf('];\ncam.ct = ['); fprintf('%f ', cam.ct);
        fprintf('];\ncam.cp = ['); fprintf('%f ', cam.cp);
        fprintf('];\n');
    end
end
