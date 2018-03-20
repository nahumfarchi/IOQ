
%paths;
% exp_dir = 'experiments/';

fp = '../../data/bunny';
[~, meshname, ~] = fileparts(fp);
mesh = MESH(fp);
X = mesh.vertices;

% [~,vid] = min( MESH.normv( [X(:,1)-min(X(:,1)) X(:,2:3)] ) );
% f = mesh.smooth_delta(vid,.4);
% % figure; MESH_VIS.func(mesh,max(f)-f,'Colormap','hot');
% 
% vf = mesh.R*(mesh.G*X(:,3));
% % figure; MESH_VIS.vf(mesh,vf);
% 
% n = 10; dt = 2/n; Dv = mesh.ifvf( vf );
% F = zeros(mesh.nv,n); F(:,1) = f;
% for i = 2:n
%     F(:,i) = expmv(-dt,Dv,F(:,i-1));
% end
% F = repmat(max(F),mesh.nv,1) - F;

% % snapshot of the simulation
% % figure; MESH_VIS.func(mesh,F(:,1),'Colormap','hot');
% % cam.pba = get(gca, 'PlotBoxAspectRatio');
% % cam.dar = get(gca, 'DataAspectRatio');
% % cam.cva = get(gca, 'CameraViewAngle');
% % cam.cuv = get(gca, 'CameraUpVector');
% % cam.ct = get(gca, 'CameraTarget');
% % cam.cp = get(gca, 'CameraPosition');
% % save([exp_dir meshname 'advection_cam.mat'],'cam');
% load([exp_dir meshname 'advection_cam.mat']);

F = [X, X];
% I = round( linspace(1,n,3) );
% MESH_IO.wfigs(meshname,mesh,F(:,I),...
%               'Colormap','hot','Montage',0,'Resolution',1024,...
%               'Camera',cam,'OpenGL',1);
MESH_IO.wfigs(meshname,mesh,F,'Montage',0,'Colormap',linspecer,'Montage',true);
          
%MESH_IO.wnrosy([meshname '.rosy'],mesh.If2v*reshape(vf,[],3));