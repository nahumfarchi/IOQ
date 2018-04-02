function create_title(filename, str, font_size)

if nargin < 3
    font_size = 14;
end

figure
%title(sprintf('$E = %.2f, ns = %d, T = %.3f$', m.miq_energy, m.n_vert_sing, elapsed));
%title(sprintf('$E = %.2f, ns = %d$', m.miq_energy, m.n_vert_sing), 'FontSize', FONT_SIZE);
title(str, 'FontSize', font_size)
axis off
set(gcf, 'WindowStyle', 'docked')
set(gcf,'color','w');

try
    export_fig(filename);
catch
    export_fig([filename '.pdf']);
end

end