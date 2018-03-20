% Plot d1*r + (pi/2)d1*p
% It might have some geometrical meaning as a function on the faces.

%% Settings
FNAME = 'sphere_s0.off';
%FNAME = 'sphere_s1.off';
%FNAME = 'ellipsoid_s4r.off';
%FNAME = 'round_cuber.off';
%FNAME = 'sphere_s0.off';
%FNAME = 'torus_fat_r2.off';
%FNAME = 'bunny.off';

%MESHES = {'sphere_s0.off', 'ellipsoid_s4r.off', 'round_cuber.off', 'torus_fat_r2.off', 'bunny.off', 'phands.off'};
MESHES = {'round_cuber.off'};

VERBOSE = true;
EPS = 1e-9;
EDGE_ALPHA = 0.2;
RESOLUTION = 1024;
AR = 1;

T = datetime;
OUT_FOLDER = fullfile('..', 'results', ...
                      [num2str(T.Year), '.', sprintf('%02d', month(T)), '.', sprintf('%02d', T.Day)], ...
                      [sprintf('%02d', T.Hour), sprintf('%02d', (T.Minute))]);
if ~exist(OUT_FOLDER, 'dir')
    mkdir(OUT_FOLDER);
end
                  
for fname = MESHES
    disp(fname)
    % Find mesh file
    p = find_data_folder();
    fp = fullfile(p, fname{:});
    
    % Load mesh
    m = Mesh();
    m.loadTM(fp);
    scale = 2 * m.avg_length;
    
    [d0, d1] = get_exterior_derivatives(m);
    N = 4;          % Direction field degree (i.e., 4 for cross field)
    f0 = [1];       % Starting face
    %v0 = [0, 0, 1]; % Starting direction
    v0 = m.V(m.F(f0, 2), :) - m.V(m.F(f0, 1), :);
    
    %% libigl MIQ (greedy rounding)
    print_header('libigl MIQ')
    res1 = nrosy_wrapper(fp, f0, v0, N);
    
    %%
    func = (d1*wrapToPi(res1.frame_diffs) + (pi/2)*d1*mod(res1.p,4));

    figure
    m.draw(func, 'EdgeAlpha', EDGE_ALPHA);
    title('d_1r+(\pi/2)d_1p')
    colorbar
    
    [~, name_only, ~] = fileparts(fname{:});
    pngname = fullfile(OUT_FOLDER, [name_only, '.png']);
    %pngname = sprintf('%s_%04d.png',filename,i);
                
    dpi = get(0, 'ScreenPixelsPerInch');
    in = RESOLUTION/dpi;

    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 AR*in in];
    fig.PaperPositionMode = 'manual';
    print(pngname, '-dpng', '-r0')
end
