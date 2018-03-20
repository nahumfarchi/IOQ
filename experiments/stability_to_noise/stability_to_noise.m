%% ========================================================================
%  Stability to noise. Add increasing noise to vertices in the normal
%  direction and test how IOQ handles it.
%  ========================================================================
FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = true;
SAVE = false;
LOAD = false;
OUT_FOLDER = 'results';
mkdir(OUT_FOLDER);

%fp = '../../../data/rounded_cube_keenan_r_3k.off';
%fp = '../../../data/bunny.off';
%fp = '../../../data/ashish_nob/dancer2.off';
fp = '../../../data/basic/cat.off';
%fp = '../../../data/ashish_nob/fandisk.off';
[~, meshname, ~] = fileparts(fp);
%OUT_FOLDER = fullfile(OUT_FOLDER, meshname);
%mkdir(OUT_FOLDER);
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;

d0 = get_exterior_derivatives(m);
L = d0'*d0;
%L = -cotmatrix(m.V, m.F);
try
    Lp = inv(gpuArray(single(full(L))) + 1/nv) - 1/nv;
catch
    Lp = invChol_mex(single(full(L)) + 1/nv) - 1/nv;
end

% Run IOQ on original mesh
[alpha1, beta1, ~, ~, out1] = IOQ_highgenus_gpu(...
                V, F, ...
                'highg_method', 'option1a', ...
                'Mesh', m, ...
                'UseGPU', USE_GPU, ...
                'InvMethod', 'ApproxResistance', ...
                'JLEps', 0.5);
            %                'LaplacianPInv', Lp, ...

k1 = [alpha1; beta1];
res1_ioq = TCODS(m, ...
    'k', k1, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);
res1_ioq.saveTM('bunny_01.ffield');
res1_ioq.saveTM('bunny_01.off');

[res1_miq, elapsed_miq] = ...
    nrosy_mex(fp, FACE0, GVEC, DEGREE);
res1_miq.saveTM('bunny_04.ffield');
res1_miq.saveTM('bunny_04.off');

results_ioq = {};
%LAMS = 0.05:0.05:0.15;
%LAMS = [0.05, 0.1, 0.15];
%LAMS = [0.05, 0.15, 0.16];
%LAMS = [0.05, 0.15, 0.25];
%LAMS = [0.05, 0.1, 0.15];
LAMS = [0.1, 0.15];
progressbar
for i = 1:length(LAMS)
    % Pertubate vertices in normal direction
    lam = LAMS(i);
    disp(lam)
    %pert = lam*m.avg_length*randn(nv, 1) .* m.VNormals;
    pert = normrnd(0, lam*m.avg_length, nv, 1);
    Vpert = V + pert;
    mp = Mesh(Vpert, F);
    

    % Run IOQ on pertubated mesh
    [alpha2, beta2, ~, ~, out2] = IOQ_highgenus_gpu(...
                    Vpert, F, ...
                    'highg_method', 'option1a', ...
                    'Mesh', mp, ...
                    'UseGPU', USE_GPU, ...
                    'InvMethod', 'ApproxResistance', ...
                    'JLEps', 0.5);
                %                    'LaplacianPInv', Lp, ...

    k2 = [alpha2; beta2];
    results_ioq{i} = TCODS(mp, ...
        'k', k2, ...
        'f0', FACE0, ...
        'theta0', THETA0, ...
        'degree', DEGREE, ...
        'CreateFField', true, ...
        'Duplicate', true, ...
        'gConstraintVec', GVEC);
    results_ioq{i}.saveTM(sprintf('bunny_0%d.ffield', i+1));
    results_ioq{i}.saveTM(sprintf('bunny_0%d.off', i+1));
    
    mp.saveTM('tmp.off');
    
    [res_miq, elapsed_miq] = nrosy_mex('tmp.off', FACE0, GVEC, DEGREE);
    results_miq{i} =  res_miq;
    results_miq{i}.saveTM(sprintf('bunny_0%d.ffield', i+4));
    results_miq{i}.saveTM(sprintf('bunny_0%d.off', i+4));
    
    progressbar(i / length(LAMS))
end

% %% Plot
% %plotting_defaults;
% close all
% cam = [];
% 
% % bunny.off
% % cam.pba = [119.5072 143.2884 119.5072];
% % cam.dar = [1 1 1];
% % cam.cva =  6.6086;
% % cam.cuv = [-0.0341 0.9948 -0.0959];
% % cam.ct  = [-0.0218 0.1097 -8.2478e-04];
% % cam.cp  = [0.1805 0.2538 1.4219];
% % POS_Y   = 0;
% 
% % cat.off
% %cam.pba = [1.4704 1.6542 1];
% %cam.dar = [1 1 1];
% %cam.cva = 6.6086;
% %cam.cuv = [0.0458 0.9771 -0.2078];
% %cam.ct  = [4.0244e-04 0.0035 4.1971e-05];
% %cam.cp  = [-3.1701 1.5100 6.3847];
% 
% %cam = [];
% %cam.pba = [1 1.1250 1.1377];
% %cam.dar = [1 1 1];
% %cam.cva = 6.6086;
% %cam.cuv = [0.0218 0.9976 -0.0656];
% %cam.ct  = [4.0244e-04 0.0035 4.1971e-05];
% %cam.cp  = [-3.8864 0.4925 6.1433];
% 
% cam = [];
% cam.pba = [1.000000 1.125000 1.137727 ];
% cam.dar = [1 1 1 ];
% cam.cva = [9.372864 ];
% cam.cuv = [0.021800 0.997600 -0.065600 ];
% cam.ct = [0.000402 0.003500 0.000042 ];
% cam.cp = [-3.886400 0.492500 6.143300 ];
% 
% POS_Y   = -0.7;
% 
% FONT_SIZE = 20;
% 
% %CM = linspecer(256);
% CM = cbrewer('seq', 'Blues', 256);
% opt = {'Func', zeros(nf, 1), 'PlotField', false', 'FaceAlpha', 1, 'EdgeAlpha', 1, 'Camera', cam, 'EdgeColor', 'none', 'Colormap', CM};
% %opt = {'PlotField', false', 'FaceAlpha', 1, 'EdgeAlpha', 1, 'EdgeColor', 'none', 'Colormap', CM};
% 
% figure
% %ha = tight_subplot(2,3,[.01 .0],0.03,0.03);
% ha = tight_subplot(2, 3, [0.01, 0.08]);
% title_ioq = sprintf('IOQ, E=%.2f, |S|=%d', res1_ioq.miq_energy, res1_ioq.n_vert_sing);
% title_miq = sprintf('MIQ, E=%.2f, |S|=%d', res1_miq.miq_energy, res1_miq.n_vert_sing);
% axes(ha(1)); res1_ioq.draw(opt{:}); th = title(title_ioq, 'FontSize', FONT_SIZE); th.Position(2) = POS_Y; material metal
% axes(ha(4)); res1_miq.draw(opt{:}); th = title(title_ioq, 'FontSize', FONT_SIZE); th.Position(2) = POS_Y; material metal
% %subplot(231); res1_ioq.draw(opt{:}); title(sprintf('IOQ, E=%.2f', res1_ioq.miq_energy))
% %subplot(234); res1_miq.draw(opt{:}); title(sprintf('MIQ, E=%.2f', res1_miq.miq_energy))
% %camlight;
% for i = 1:2
%     title_ioq = sprintf('E=%.2f, |S|=%d', results_ioq{i}.miq_energy, results_ioq{i}.n_vert_sing);
%     title_miq = sprintf('E=%.2f, |S|=%d', results_miq{i}.miq_energy, results_miq{i}.n_vert_sing);
%     
%     axes(ha(i+1)); results_ioq{i}.draw(opt{:}); material metal
%     %th = title(sprintf('\\sigma(0, %.2e), E=%.2f, |S|=%d', LAMS(i)*m.avg_length, results_ioq{i}.miq_energy, results_ioq{i}.n_vert_sing));
%     th = title(title_ioq, 'FontSize', FONT_SIZE);
%     th.Position(2) = POS_Y;
%     axes(ha(i+4)); results_miq{i}.draw(opt{:}); material metal
%     %th = title(sprintf('\\sigma(0, %.2e), E=%.2f, |S|=%d', LAMS(i)*m.avg_length, results_miq{i}.miq_energy, results_miq{i}.n_vert_sing))
%     th = title(title_miq, 'FontSize', FONT_SIZE);
%     th.Position(2) = POS_Y;
%     
%     %subplot(2,3,i+1); results_ioq{i}.draw(opt{:})
%     %title(sprintf('sigma(0, %.3g), E=%.2f', LAMS(i)*m.avg_length, results_ioq{i}.miq_energy))
%     %subplot(2,3,i+4); results_miq{i}.draw(opt{:})
%     %title(sprintf('sigma(0, %.3g), E=%.2f', LAMS(i)*m.avg_length, results_miq{i}.miq_energy))
%     %camlight
% end
% 
% set(gcf, 'WindowStyle', 'docked')
% set(gcf,'color','w');
% 
% export_fig('stability_to_noise.pdf')
% export_fig('stability_to_noise.png')
% 
% if SAVE
%     filename = fullfile(OUT_FOLDER, [meshname '_stability_to_noise']);
%     print(gcf, filename, '-dpng', '-r300')
%     saveas(gcf, filename, 'fig')
% end


%% Plot
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaulttextinterpreter','latex')

close all
cam = [];

% bunny.off
% cam.pba = [119.5072 143.2884 119.5072];
% cam.dar = [1 1 1];
% cam.cva =  6.6086;
% cam.cuv = [-0.0341 0.9948 -0.0959];
% cam.ct  = [-0.0218 0.1097 -8.2478e-04];
% cam.cp  = [0.1805 0.2538 1.4219];
% POS_Y   = 0;

% cat.off
%cam.pba = [1.4704 1.6542 1];
%cam.dar = [1 1 1];
%cam.cva = 6.6086;
%cam.cuv = [0.0458 0.9771 -0.2078];
%cam.ct  = [4.0244e-04 0.0035 4.1971e-05];
%cam.cp  = [-3.1701 1.5100 6.3847];

%cam = [];
%cam.pba = [1 1.1250 1.1377];
%cam.dar = [1 1 1];
%cam.cva = 6.6086;
%cam.cuv = [0.0218 0.9976 -0.0656];
%cam.ct  = [4.0244e-04 0.0035 4.1971e-05];
%cam.cp  = [-3.8864 0.4925 6.1433];

cam = [];
cam.pba = [1.000000 1.125000 1.137727 ];
cam.dar = [1 1 1 ];
cam.cva = [9.372864 ];
cam.cuv = [0.021800 0.997600 -0.065600 ];
cam.ct = [0.000402 0.003500 0.000042 ];
cam.cp = [-3.886400 0.492500 6.143300 ];

POS_Y   = -0.6;

FONT_SIZE = 20;

%CM = linspecer(256);
%CM = cbrewer('seq', 'Blues', 256);
CM = cbrewer('div', 'RdBu', 256);
%opt = {'Func', zeros(nf, 1), 'PlotField', false', 'FaceAlpha', 1, 'EdgeAlpha', 1, 'Camera', cam, 'EdgeColor', 'none', 'Colormap', CM};
%opt = {'PlotField', false', 'FaceAlpha', 1, 'EdgeAlpha', 1, 'EdgeColor', 'none', 'Colormap', CM};
opt = {'PlotField', false, ...
       'PlotSing', true, ...
       'FaceAlpha', 1, ...
       'EdgeAlpha', 0, ...
       'Dock', true, ...
       'View', VIEW, ...
       'Colormap', CM, ...
       'Caxis', 'auto', ...
       'Camera', cam};

%figure
%ha = tight_subplot(2,3,[.01 .0],0.03,0.03);
%ha = tight_subplot(2, 3, [0.01, 0.08]);
title_ioq = sprintf('IOQe, $E=%.2f, |S|=%d$', res1_ioq.miq_energy, res1_ioq.n_vert_sing);
title_miq = sprintf('MIQ, $E=%.2f, |S|=%d$', res1_miq.miq_energy, res1_miq.n_vert_sing);

figure; res1_ioq.draw(opt{:}); th = title(title_ioq, 'FontSize', FONT_SIZE); th.Position(2) = POS_Y; material shiny
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
export_fig('stability_to_noise_01.pdf');

figure; res1_miq.draw(opt{:}); th = title(title_miq, 'FontSize', FONT_SIZE); th.Position(2) = POS_Y; material shiny
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');
export_fig('stability_to_noise_04.pdf');

%subplot(231); res1_ioq.draw(opt{:}); title(sprintf('IOQ, E=%.2f', res1_ioq.miq_energy))
%subplot(234); res1_miq.draw(opt{:}); title(sprintf('MIQ, E=%.2f', res1_miq.miq_energy))
%camlight;
for i = 1:2
    title_ioq = sprintf('$E=%.2f, |S|=%d$', results_ioq{i}.miq_energy, results_ioq{i}.n_vert_sing);
    title_miq = sprintf('$E=%.2f, |S|=%d$', results_miq{i}.miq_energy, results_miq{i}.n_vert_sing);
    
    figure; results_ioq{i}.draw(opt{:}); material shiny %metal
    %th = title(sprintf('\\sigma(0, %.2e), E=%.2f, |S|=%d', LAMS(i)*m.avg_length, results_ioq{i}.miq_energy, results_ioq{i}.n_vert_sing));
    th = title(title_ioq, 'FontSize', FONT_SIZE);
    th.Position(2) = POS_Y;
    set(gcf, 'WindowStyle', 'docked')
    set(gcf,'color','w');
    export_fig(sprintf('stability_to_noise_0%d.pdf', i+1));
    
    figure; results_miq{i}.draw(opt{:}); material shiny %metal
    %th = title(sprintf('\\sigma(0, %.2e), E=%.2f, |S|=%d', LAMS(i)*m.avg_length, results_miq{i}.miq_energy, results_miq{i}.n_vert_sing))
    th = title(title_miq, 'FontSize', FONT_SIZE);
    th.Position(2) = POS_Y;
    set(gcf, 'WindowStyle', 'docked')
    set(gcf,'color','w');
    export_fig(sprintf('stability_to_noise_0%d.pdf', i+4));
    
    %subplot(2,3,i+1); results_ioq{i}.draw(opt{:})
    %title(sprintf('sigma(0, %.3g), E=%.2f', LAMS(i)*m.avg_length, results_ioq{i}.miq_energy))
    %subplot(2,3,i+4); results_miq{i}.draw(opt{:})
    %title(sprintf('sigma(0, %.3g), E=%.2f', LAMS(i)*m.avg_length, results_miq{i}.miq_energy))
    %camlight
end

set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');

%export_fig('stability_to_noise.pdf')
%export_fig('stability_to_noise.png')

% if SAVE
%     filename = fullfile(OUT_FOLDER, [meshname '_stability_to_noise']);
%     print(gcf, filename, '-dpng', '-r300')
%     saveas(gcf, filename, 'fig')
% end



set(groot, 'defaultAxesTickLabelInterpreter','none'); 
set(groot, 'defaultLegendInterpreter','none');
set(0,'defaulttextinterpreter','none')