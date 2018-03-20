NAME = 'squares4_row.off';
N = 4;
p0 = [0, 0, 0];
p1 = [0, 1, 0];
p0_vid = 1;
p1_vid = 2;
next_vid = 3;
F = [];
V = [p0; p1];
% p1 ---- q1
% | *    * |
% |  mid   |
% | *    * |
% p0 ---- q0
for i = 1:N
    q0 = [p0(1)+1, p0(2), 0];
    q1 = [p1(1)+1, p1(2), 0];
    mid = [p0(1)+0.5, p0(2)+0.5, 0];
    
    q1_vid = next_vid;
    q0_vid = next_vid+1;
    mid_vid = next_vid+2;
    next_vid = next_vid+3;
    
    f1 = [p0_vid, mid_vid, p1_vid];
    f2 = [p1_vid, mid_vid, q1_vid];
    f3 = [q1_vid, mid_vid, q0_vid];
    f4 = [q0_vid, mid_vid, p0_vid];
    
    V = [V; q1; q0; mid];
    F = [F; f1; f2; f3; f4];
    
    p0 = q0;
    p1 = q1;
    p0_vid = q0_vid;
    p1_vid = q1_vid;
end

m = Mesh(V, F);
figure;
m.draw();
hold on;
m.drawLabels(true, true, true);
hold off;
view([0,0,1])
axis off
p = find_data_folder();
file_path = fullfile(p, NAME);
m.saveTM(file_path);