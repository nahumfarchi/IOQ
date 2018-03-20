function plot_resistance(m, R, v, varargin)
%plot_resistance(m, R)

if nargin < 3
    v = 1;
end

try
    hold on
    f = R(:, v);
    colors = linspecer(2);
    H = plot3( m.V(v,1), m.V(v,2), m.V(v,3), '.', 'markers', 40, 'color', colors(2,:));
    m.draw(f, 'FaceAlpha', 1, 'EdgeAlpha', 0, varargin{:});
    hold off
catch E
    hold off
    throw(E);
end

