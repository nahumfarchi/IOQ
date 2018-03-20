B = [2, -1, 0, -1; -1, 2, 0, -1; 0, 0, 2, -2];
Q = B'*B;

% Shortest vector
G = -Q;
G = G - diag(diag(G));

[mv, mw] = min_cut(G);

v1 = sum(B(:,mv==1), 2);
v2 = sum(B(:,mv==0), 2);

uG = graph(G);
[mf,~,cs,ct] = maxflow(uG,1,size(G,1));

% Closest point
x1 = [-4; 2; 0];
x2 = [-5; 1; 0];
y = x1 + (x2-x1)/2;

y(1) = y(1)+0.1;
x = closest_lpoint_vor(B, y);
check_norm('x', 'x1');

y(2) = y(2)-0.2;
x = closest_lpoint_vor(B, y);
check_norm('x', 'x2');

y = [-3.1; 0; 0];
x = closest_lpoint_vor(B, y);
check_norm('x', '[-3;0;0]');
