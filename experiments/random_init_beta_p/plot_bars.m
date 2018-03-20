function h = plot_bars(datas, lgds, xlabels)
%PLOT_BARS Summary of this function goes here
%   Detailed explanation goes here
    y_err = [];
    y = [];
    for i = 1:numel(datas)
       data = datas{i};
       [y_err1, y1] = get_bar_data(data);
       y_err = [y_err; y_err1];
       y = [y; y1];
    end

    h = barwitherr(y_err, y);
    legend(lgds);
    ylim([-10, max(max(y))+20])
    set(gca, 'XTickLabel', xlabels)
    for i = 1:length(h(1).XData)
        text(h(1).XData(i), h(1).YData(i), ['Emiq=', num2str(h(1).YData(i))], ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom')
    end
    
% subplot(224)
%     [y_err1, y1] = get_bar_data(datas{1});
%     [y_err2, y2] = get_bar_data(datas{2});
%     [y_err3, y3] = get_bar_data(datas{3});
%     y_err = [y_err1; y_err2; y_err3];
%     y = [y1; y2; y3];    
%     
%     h = barwitherr(y_err, y);
%     legend('Emiq', 'E1', 'E2', 'E', 'm')
%     set(gca, 'XTickLabel', {'Opt 3 rand beta', 'Opt 3 beta=round', 'opt 3 beta=0'})
%     for i = 1:length(h(1).XData)
%         text(h(1).XData(i), h(1).YData(i), ['Emiq=', num2str(h(1).YData(i))], ...
%             'HorizontalAlignment', 'center', ...
%             'VerticalAlignment', 'bottom')
%     end
%     title('Option 3 final energy')
