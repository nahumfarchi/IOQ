function [res] = panorama(imgs, out, direction, irfanview_path)
    %function [res] = panorama(imgs, out, direction, irfanview_path)
    %
    %Create a panorama image from the given images.
    % Input:
    %   imgs - a list of images
    %   out  - output filename
    %   direction - pan direction ('horizontal' or 'vertical')
    %   irfanview_path - path to irfanview exe (optional)
    
    if nargin < 3
        direction = "horizontal";
    end
    
    switch direction
        case "horizontal"
            direction = 1;
        case "vertical"
            direction = 2;
    end
    
    if nargin < 4
        irfanview_path = "C:/Program Files (x86)/IrfanView/i_view32.exe";
    end

    cmd = sprintf('"%s" /panorama=(%d', irfanview_path, direction);
    
    for j=1:length(imgs)
        cmd = [cmd ',' imgs{j}];
    end
    cmd = [cmd ') /dpi=(300,300) /spacing=100 /convert=' out];
    res = system(cmd);
end

