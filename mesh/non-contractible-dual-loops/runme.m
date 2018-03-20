clear all; close all

addpath(genpath('.'));

meshname = 'eight';
mesh = MESH( meshname );
nf = mesh.nf; nv = mesh.nv;

%%%%% Primal/Dual tree-co-tree decomposition. Assumes no boundaries!

% Primal graph
s = reshape(mesh.triangles,[],1);
t = reshape(mesh.triangles(:,[2,3,1]),[],1);
Gp = digraph(s,t); 
Ap = adjacency(Gp); 
Gp = graph(Ap);

ne = sum(sum(Ap));
[ii,jj] = find(Ap);
Ep = sparse(ii,jj,1:ne,nv,nv);

% Dual graph
F1 = sparse(s,t,[1:nf,1:nf,1:nf]',nv,nv);
F2 = sparse(t,s,[1:nf,1:nf,1:nf]',nv,nv);
[ii,jj,ss] = find(Ep); ll = sub2ind(size(Ep),ii,jj);
Ed = sparse(F1(ll),F2(ll),ss,nf,nf);
Ad = double(Ed ~= 0);
Gd = graph(Ad);

if (sum(sum(Ad)) ~= sum(sum(Ap)))
    error('?');
end

% Primal spanning tree
Tp = minspantree(Gp,'Method','sparse'); Apt = adjacency(Tp);
% Primal edges not in tree
Apnt = double(Ap & ~Apt);
[ii,jj] = find(Apnt); ll = sub2ind(size(Ep),ii,jj);
%figure; MESH_VIS.mesh(mesh); hold on; 
%[ii,jj] = find(Apt); line([mesh.vertices(ii,1),mesh.vertices(jj,1)]',[mesh.vertices(ii,2),mesh.vertices(jj,2)]',[mesh.vertices(ii,3),mesh.vertices(jj,3)]','linewidth',2,'color','b'); 

% Dual graph without primal tree
Ednpt = sparse(F1(ll),F2(ll),Ep(ll),nf,nf);
Adnpt = double(Ednpt ~= 0);
Gdnpt = graph(Adnpt);

% dual spanning tree
[Td, pred] = minspantree(Gdnpt,'method','sparse'); Adt = adjacency(Td);
Xf = mesh.Iv2f*mesh.vertices;
%[ii,jj] = find(Adt); line([Xf(ii,1),Xf(jj,1)]',[Xf(ii,2),Xf(jj,2)]',[Xf(ii,3),Xf(jj,3)]','linewidth',2,'color','r'); 

% num edges in dual graph not in primal spanning tree
nednpt = sum(sum(Adnpt));
% num edges in dual spanning tree
nedt = sum(sum(Adt));
% number of edges which are in neither trees = 2*g
if ((nednpt - nedt)/2 ~= -(nv + nf - ne/2 - 2))
    error('?');
end

% Edges in neither trees
R = Adnpt - Adt;
if (norm(R - R','fro')~=0)
    error('?');
end
R = R - triu(R);
[f1s,f2s] = find(R);
%[ii,jj] = find(R); line([Xf(ii,1),Xf(jj,1)]',[Xf(ii,2),Xf(jj,2)]',[Xf(ii,3),Xf(jj,3)]','linewidth',4,'color','m'); 

cycles = {};
for i=1:length(f1s)
    figure; MESH_VIS.mesh(mesh); hold on; 
    line([Xf(f1s(i),1),Xf(f2s(i),1)]',[Xf(f1s(i),2),Xf(f2s(i),2)]',[Xf(f1s(i),3),Xf(f2s(i),3)]','linewidth',4,'color','m'); 
    cycle_a = f1s(i);
    while (pred(cycle_a(end))~=0)
        cycle_a = [cycle_a, pred(cycle_a(end))];
        v1 = cycle_a(end-1); v2 = cycle_a(end);
        line([Xf(v1,1),Xf(v2,1)]',[Xf(v1,2),Xf(v2,2)]',[Xf(v1,3),Xf(v2,3)]','linewidth',4,'color','b'); 
    end
    cycle_b = f2s(i);
    while (pred(cycle_b(end))~=0)
        cycle_b = [cycle_b, pred(cycle_b(end))];
        v1 = cycle_b(end-1); v2 = cycle_b(end);
        line([Xf(v1,1),Xf(v2,1)]',[Xf(v1,2),Xf(v2,2)]',[Xf(v1,3),Xf(v2,3)]','linewidth',4,'color','r'); 
    end
    cycles{end+1} = [cycle_a(end:-1:1),cycle_b];
end




