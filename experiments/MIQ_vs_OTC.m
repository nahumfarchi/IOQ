
if ~exist('experiment_name', 'var')
    disp('Set variable experiment_name!')
    return
end

print_header(['experiment_name : ', experiment_name])

%experiment_name = 'bunnies_connectivity_lap';
base_folder = '../results/experiments';
data_folder = fullfile(base_folder, experiment_name, 'ffields');
out_folder = fullfile(base_folder, experiment_name, 'field_plots_N4');
if ~exist(base_folder, 'dir')
    error('Base folder does not exist')
end
if ~exist(data_folder, 'dir')
    error('Data folder does not exist')
end
if ~exist(out_folder, 'dir')
    mkdir(out_folder);
end

overwrite = false;

resolution = 4048;

filepaths = get_filepaths(data_folder, '.ffield', false);
i = 1;
while i < numel(filepaths)
%for i = 1 : 2 : numel(filepaths) - 1
    fp_miq = filepaths{i};
    fp_otc = filepaths{i+1};
    
    [~, name1, ~] = fileparts(fp_miq);
    [~, name2, ~] = fileparts(fp_otc);
    name1 = strsplit(name1, '_');
    name1 = name1{1};
    name2 = strsplit(name2, '_');
    name2 = name2{1};
    
    if ~strcmp(name1, name2)
        disp([name1, ' and ', name2, ' do not match. Skipping...'])
        i = i + 1;
        continue
    end
    
    assert(strfind(fp_miq, 'MIQ.ffield'))
    assert(strfind(fp_otc, 'OTC_blockinv_gpuloop.ffield'))
    
    out_path = fullfile(out_folder, [name1, '.png']);
    plot_exists = exist(out_path, 'file');
    if ~overwrite && plot_exists 
        disp(['Skipping ', name1, ' (already exists)...'])
        i = i + 2;
        continue
    elseif overwrite && plot_exists
        disp(['Overwriting ', out_path, '...'])
    end
    
    disp(['Plotting ', name1, '...'])
    
    m_miq = Mesh();
    m_otc = Mesh();
    m_miq.loadTM(fp_miq);
    m_otc.loadTM(fp_otc);
    
    fig = figure();
    ha = tight_subplot(3,2,[.01 0],[.1 .1],[0 0]);
   %for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
    
    
    %subplot(321)
    axes(ha(1));
    m_miq.draw()
    title({['miq, E = ', num2str(m_miq.miq_energy)], ['singularities: ', num2str(m_miq.n_singularities)]})
    view(2);
    
    %subplot(322)
    axes(ha(2));
    m_otc.draw()
    title({['otc, E = ', num2str(m_otc.miq_energy)], ['singularities: ', num2str(m_otc.n_singularities)]})
    view(2)
    
    %subplot(323)
    axes(ha(3));
    m_miq.draw()
    view(3)
    
    %subplot(324)
    axes(ha(4));
    m_otc.draw()
    view(3)
    
    %subplot(325)
    axes(ha(5));
    m_miq.draw()
    view(180, -90)
    
    %subplot(326)
    axes(ha(6));
    m_otc.draw()
    view(180, -90)
    
    %set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
    
    mysave(fig, name1, out_folder, resolution)
    %export_fig('test.png')
    
    i = i + 2;
end 