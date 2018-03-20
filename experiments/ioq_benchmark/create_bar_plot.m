clear all
close all

stats_file = 'results\ashish_nob\stats_clean';

s = load(stats_file); s = s.stats_clean;

for i=2:size(s,1)-2
    names{i-1} = s{i,1};
    for j=2:size(s,2)
        d = s{i,j};
        T(i-1,j-1) = d.elapsed_total;
        E(i-1,j-1) = d.miq_energy;
        S(i-1,j-1) = d.n_singularities;
    end
end

[ss,is] = sort(E(:,2)+T(:,2));
S = S(is,:);
T = T(is,:);
E = E(is,:);
names = names(is);

for i=2:3
    impr(:,i-1) = E(:,1).*(S(:,1)+1)./E(:,i)./(S(:,i)+1);
    imprT(:,i-1) = T(:,1)./T(:,i);
end

col = linspecer(3); cc = col(2,:); col(2,:) = col(3,:); col(3,:) = cc;

figure; hold on;
c = categorical(names);
h = barh(c,impr); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
legend('(E*ns)_M / (E*ns)_I','(E*ns)_M / (E*ns)_J');
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
xlim([0 7]);
set(gcf, 'WindowStyle', 'docked')

figure; hold on;
c = categorical(names);
h = barh(c,imprT); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
legend('T_M / T_I','T_M / T_J');
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
xlim([0 1.2]);
set(gcf, 'WindowStyle', 'docked')

figure; hold on;
c = categorical(names);
h = barh(c,S(:,[2,3,1])); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
h(3).FaceColor = col(3,:);
legend('#s IOQ', '#s IOQ_J','#s MIQ');
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
xlim([0 50]);
set(gcf, 'WindowStyle', 'docked')

figure; hold on;
c = categorical(names);
h = barh(c,E(:,[2,3,1])); set(gca, 'TickLabelInterpreter', 'none');
h(1).FaceColor = col(1,:);
h(2).FaceColor = col(2,:);
h(3).FaceColor = col(3,:);
legend('E IOQ', 'E IOQ_J','E MIQ');
ax = gca; set(gca,'XGrid','on','Ydir','reverse');
xlim([0 500]);
set(gcf, 'WindowStyle', 'docked')

