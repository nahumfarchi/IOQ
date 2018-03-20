clear all
close all

stats_file = 'results\ashish_nob\stats_clean';

s = load(stats_file); s = s.stats_clean;

for i=2:size(s,1)-2
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

nm = 20;
col = linspecer(nm);
models = 1:nm; sf = 250;
figure; hold on;

% A triangle per mesh
for i=1:length(models)
    x = [S(models(i),:); E(models(i),:)]';
    patch('faces',[1,2,3],'vertices',x,'facecolor',col(i,:),'edgecolor','k','facealpha',.3);
end

xlabel('sing');
ylabel('E');

% Square has a smaller marker than circle and diamond
scatter(S(models,2),E(models,2),T(models,2)*sf,col,'diamond','filled','MarkerFaceAlpha',.7,'MarkerEdgeColor','k'); 
scatter(S(models,3),E(models,3),T(models,3)*sf*7/5,col,'square','filled','MarkerFaceAlpha',.7,'MarkerEdgeColor','k'); 
scatter(S(models,1),E(models,1),T(models,1)*sf,col,'filled','MarkerFaceAlpha',.7,'MarkerEdgeColor','k'); 

return

% % A triangle per mesh
% models1 = models(1:2:end); models2 = models(2:2:end); 
% col1 = col(1:2:end,:);col2 = col(2:2:end,:);
% for i=1:length(models1)
%     x = [S(models1(i),:); E(models1(i),:)]';
%     patch('faces',[1,2,3],'vertices',x,'facecolor',col1(i,:),'edgecolor','k','facealpha',.3);
% end
% for i=1:length(models2)
%     x = [S(models2(i),:); E(models2(i),:)]';
%     patch('faces',[1,2,3],'vertices',x,'facecolor',col2(i,:),'edgecolor','k','facealpha',.3);
% end
% 
% % Square has a smaller marker than circle and diamond
% scatter(S(models2,2),E(models2,2),T(models2,2)*sf,col2,'diamond','filled','MarkerFaceAlpha',.7,'MarkerEdgeColor','k'); 
% scatter(S(models1,2),E(models1,2),T(models1,2)*sf,col1,'diamond','filled','MarkerFaceAlpha',.7,'MarkerEdgeColor','k'); 
% scatter(S(models2,3),E(models2,3),T(models2,3)*sf*7/5,col2,'square','filled','MarkerFaceAlpha',.7,'MarkerEdgeColor','k'); 
% scatter(S(models1,3),E(models1,3),T(models1,3)*sf*7/5,col1,'square','filled','MarkerFaceAlpha',.7,'MarkerEdgeColor','k'); 
% scatter(S(models2,1),E(models2,1),T(models2,1)*sf,col2,'filled','MarkerFaceAlpha',.7,'MarkerEdgeColor','k'); 
% scatter(S(models1,1),E(models1,1),T(models1,1)*sf,col1,'filled','MarkerFaceAlpha',.7,'MarkerEdgeColor','k'); 

% scatter(S(models,2)+1,E(models,2),T(models,2)*sf+1,col,'diamond','filled','MarkerEdgeColor','k'); 
% scatter(S(models1,2)+1,E(models1,2),T(models1,2)*(sf+1)*.7,col1,'diamond','filled','MarkerEdgeColor','k'); 
% scatter(S(models,3)+1,E(models,3),(T(models,3)*sf+1)*7/5,col,'square','filled','MarkerEdgeColor','k'); 
% scatter(S(models1,3)+1,E(models1,3),(T(models1,3)*sf+1)*7/5*.7,col1,'square','filled','MarkerEdgeColor','k'); 
% scatter(S(models,1)+1,E(models,1),T(models,1)*sf+1,col,'filled','MarkerEdgeColor','k'); 
% scatter(S(models1,1)+1,E(models1,1),T(models1,1)*(sf+1)*.7,col1,'filled','MarkerEdgeColor','k'); 
% 
% scatter(S(models,2)+1,E(models,2),30,col,'diamond','filled','MarkerEdgeColor','k'); 
% scatter(S(models1,2)+1,E(models1,2),30*.7,col1,'diamond','filled','MarkerEdgeColor','k'); 
% scatter(S(models,3)+1,E(models,3),30*7/5,col,'square','filled','MarkerEdgeColor','k'); 
% scatter(S(models1,3)+1,E(models1,3),30*7/5*.7,col1,'square','filled','MarkerEdgeColor','k'); 
% scatter(S(models,1)+1,E(models,1),30,col,'filled','MarkerEdgeColor','k'); 
% scatter(S(models1,1)+1,E(models1,1),30*.7,col1,'filled','MarkerEdgeColor','k'); 
% scatter(S(1:nm,1)+1,E(1:nm,1),50,col,'filled'); 
% scatter(S(1:nm,2)+1,E(1:nm,2),50,col,'diamond','filled'); 
% scatter(S(1:nm,3)+1,E(1:nm,3),70,col,'square','filled'); 
