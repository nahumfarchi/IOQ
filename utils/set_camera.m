function set_camera(ca,cam)
    if isempty(cam)
        view(3);
        return
    end
    set(ca, 'PlotBoxAspectRatio',cam.pba);
    set(ca, 'DataAspectRatio',cam.dar);
    set(ca, 'CameraViewAngle',cam.cva);
    set(ca, 'CameraUpVector',cam.cuv);
    set(ca, 'CameraTarget',cam.ct);
    set(ca, 'CameraPosition',cam.cp);
end
