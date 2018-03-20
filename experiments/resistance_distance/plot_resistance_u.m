function plot_resistance_u(m, R, u, vert)
%plot_resistance(m, R)

if nargin < 4
    vert = 1;
end

try
    hold on
    f = R(:, vert) - u(vert) + u;
    %f = -u(vert) + u;
    colors = linspecer(2);
    H = plot3( m.V(vert,1), m.V(vert,2), m.V(vert,3), '.', 'color', colors(2,:));
    m.draw(f, 'FaceAlpha', 1, 'EdgeAlpha', 0);
    hold off
catch E
    hold off
    throw(E);
end

