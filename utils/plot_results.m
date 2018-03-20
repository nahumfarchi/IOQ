function plot_results(plots, results, legends, out_folder)
    % Input:
    %   plots - cellarray of structures, each describing what to plot. Each
    %   struct has the following fields:
    %       title
    %       name
    %       get_x_field
    %       get_y_field
    %       style
    %   results - cellarray of results returned by run_experiments
    %   legends
    %   out_folder
    
    is_valid = @(x) ~isempty(x) && ...
                    (~isfield(x, 'failed') || ...
                    ~x.failed);
    legends = {};
    for i = 1:numel(plots)
        
        p = plots{i};
        hold on
        
        if isfield(p, 'y_cols')
            xx = cellfun(p.get_x_field, ...
                results(:, p.x_col), ...
                'UniformOutput', false);
            [xx, sort_inds] = sort(cell2mat(xx));
            
            data = results(:, p.y_cols);
            valid = cellfun(@(x) is_valid(x), data);
            [rows, cols] = size(data);
            yy = inf(rows, cols);
            yy(valid) = cellfun(p.get_y_field, data(valid));
            [yy, col_inds] = min(yy, [], 2);
            
            yy = yy(sort_inds);
            col_inds = col_inds(sort_inds);
            
            plot(xx, yy, p.style{:})
            legends{end+1} = p.legend;
            
            for c = p.y_cols
                %c = col_inds(i);
                idx = col_inds == c;
                if nnz(idx) > 0
                    scatter(xx(idx), yy(idx), p.symbols{c});
                    legends{end+1} = p.symbol_legends{c};
                end
            end
            
            
        else
            figs(i) = figure();
            cols = size(results, 2);
            
            for j = 1:cols
                yy = results(:, j);
                valid = cellfun(@(x) isvalid(x), yy);
                %suc = cellfun(@(x) ~x.failed, yy(valid));
                yy = yy(valid);
                xx = cellfun(p.get_x_field, yy);
                yy = cellfun(p.get_y_field, yy);

                [xx, inds] = sort(xx);
                yy = yy(inds);

                plot(xx, yy, p.style{:})
            end

            legend(legends{:}, 'Location', 'northwest')
            title(p.title)
            mysave(figs(i), p.name, out_folder);
        end
        
        hold off
        
    end
    
    legend(legends)
    
end