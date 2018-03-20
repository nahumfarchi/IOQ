function mysave(fig, filename, out_folder, resolution, AR)
    % function mysave(fig, filename, out_folder)
    if nargin < 4
        resolution = 1024;
    end
    if nargin < 5
        AR = 1;
    end
    
    %pngname = sprintf('%s_%04d.png', filename, i);
    filename = fullfile(out_folder, filename);
    pngname = sprintf('%s.png', filename);

    dpi = get(0, 'ScreenPixelsPerInch');
    in = resolution/dpi;

    %fig = gcf;
    figure(fig);
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 AR*in in];
    fig.PaperPositionMode = 'manual';
    print(pngname,'-dpng','-r0')
end