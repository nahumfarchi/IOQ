%%
filepaths = get_filepaths(OUT_FOLDER_LP, '.mat');
n_files = numel(filepaths);
sparsity = [];
n_verts = [];
T_load = [];
T_inv = [];

for i = 1:n_files
    fp = filepaths{i};
    print_header(sprintf('%d // %d', i, n_files))
    disp(fp)
    
    tic
    load(fp);
    T_load(end+1) = toc;
    
    T_inv(end+1) = elapsed_inv;
    
    n_verts(end+1) = size(Lp, 1);
    sparsity(end+1) = length(find(Lp>1e-10)) / size(Lp,1)^2;
    
end

%% Sort results by vertex size
[n_verts, I] = sort(n_verts);
sparsity = sparsity(I);
T_load = T_load(I);
T_inv = T_inv(I);

% Plotting defaults
width = 3;     % Width in inches
height = 3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

% The new defaults will not take effect if there are any open figures. To
% use them, we close all figures, and then repeat the first example.
close all;

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);

figure(1);
plot(n_verts, sparsity, 'X');
xlabel('# verts')
ylabel('sparsity %')
title('Sparsity % of L^{\dagger}')
print('Lp_sparsity', '-dpng', '-r300');

figure(2); 
plot(n_verts, T_load, 'X', n_verts, T_inv, 'X')
xlabel('# verts')
ylabel('Time')
legend('Load time', 'Inv time')
title('Time to load Lp from disk vs time to invert L')
print('Lp_time', '-dpng', '-r300');

