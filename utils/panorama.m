function [res] = panorama(imgs, out, direction, spacing, irfanview_path)
    %function [res] = panorama(imgs, out, direction, spacing, irfanview_path)
    %
    %Create a panorama image from the given images.
    % Input:
    %   imgs - a list of images
    %   out  - output filename
    %   direction - pan direction ('horizontal' or 'vertical')
    %   spacing - spacing between images in ems (default is 2em). TODO
    %   irfanview_path - path to irfanview exe (optional)
    
    em = 16; % 1em is 16px
    if nargin < 3
        direction = "horizontal";
    end
%     if nargin < 4
%         spacing = 1;
%     end

    switch direction
        case "horizontal"
            direction = 1;
        case "vertical"
            direction = 2;
    end
    
    if nargin < 5
        irfanview_path = "C:/Program Files (x86)/IrfanView/i_view32.exe";
    end
    
%     if spacing > 0
%         blank_image = ones(spacing, spacing, 3);
%         blank_name = [tempname, '.png'];
%         imwrite(blank_image, blank_name);
%     end

    cmd = sprintf('"%s" /panorama=(%d', irfanview_path, direction);    
    
    for j=1:length(imgs)
        cmd = [cmd ',' imgs{j}];
%         if spacing > 0 && j ~= length(imgs)
%             cmd = [cmd ',' blank_name];
%         end
    end
    
    cmd = [cmd ') /dpi=(300,300) /convert=' out];
    res = system(cmd);
    
%     if spacing > 0
%         delete(blank_name);
%     end
end

