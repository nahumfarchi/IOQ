%to run this script without crashing matlab, use:
%cmd_str='matlab -nosplash -nodesktop -r "setup; conformal_seamless; exit" &';
%system(cmd_str)
%or:
%cmd_str='matlab -nosplash -nodesktop -r "setup; waitforbuttonpress; conformal_seamless; exit" &';
%21.1.17: not relevant anymore since it stopped crashing after
% updating OpenMesh.
%%
%experiments_setup;

meshname = '3holes';
%meshname = 'torus_fat_r2';
out_folder = './results/';
inext = '.off';
outext = '.obj';
if ~exist(out_folder, 'dir')
    mkdir(out_folder);
end
fp_in = fullfile('..', '..', '..', 'data', [meshname inext]);
uv_scale_factor = 10;

m = Mesh(fp_in);
V = m.V;
F = m.F;
nv = m.nV; nf = m.nF;
cycles = m.generator_cycles();
ng2 = numel(cycles);
[A, K, d0, d1, H] = tcods_gsystem(V, F);

kioq = IOQ_cpu(V, F);
alpha_p_ioq = (4 - kioq) * (pi / 2); % alpha_hat in marcel's code
beta_p_ioq = round(-K(end-1:end) / (pi/2)) * pi/2; % kappa_hat in marcel's code

[theta, p, kmiq, R, local_frames, frame_diffs, elapsed, pFixed] = ...
    NRoSy_mex_bin(fp_in, [1], [1,0,0], 4);
alpha_p_miq = (4 - (-4*kmiq)) * (pi / 2);
beta_p_miq = ((K(end-1:end) - H'*frame_diffs) / (2*pi) - 0.25*H'*p) * (4*pi/2);
%tmp = (K(end-1:end) - H'*frame_diffs) / (2*pi) + 0.25*H'*p;
%beta_p_miq = (4 - (-4*


%inds = randperm(nv, 8);
%theta_hat = 2*pi*ones(nv, 1);
%theta_hat(inds) = 3*pi/2;


%kappa_hat = [];

%Create a vector of cycles (face-id based), separated by -1
gamma = [];
if ng2 > 0
    for i = 1:ng2-1
        cy = cycles{i};
        gamma = [gamma; cy(:); -1];
    end
    cy = cycles{end};
    gamma = [gamma; cy(:)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 22.11.17
%% IOQ
fp_out = fullfile(out_folder, [meshname '_ioq' outext]);
[u, v] = conformal_seamless_mex_bin(...
    fp_in, fp_out, alpha_p_ioq, [], gamma, uv_scale_factor);

muv = read_obj(fp_out);
V = [muv.vt, zeros(size(muv.vt, 1), 1)]; F = muv.f.vt;
write_obj(out_folder, [meshname '_ioq_uv' outext], V, F);

%% MIQ
fp_out = fullfile(out_folder, [meshname '_miq' outext]);
[u, v] = conformal_seamless_mex_bin(...
    fp_in, fp_out, alpha_p_miq, beta_p_miq, gamma, uv_scale_factor);

muv = read_obj(fp_out);
V = [muv.vt, zeros(size(muv.vt, 1), 1)]; F = muv.f.vt;
write_obj(out_folder, [meshname '_miq_uv' outext], V, F);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 21.11.17
% %% IOQ
% theta_hat = (4 - kioq) * (pi / 2);
% kappa_hat = round(flatten_generators(m) / (pi/2)) * pi/2;
% fp_out = fullfile(out_folder, [meshname '_ourcycles_ioq' outext]);
% [u, v] = conformal_seamless_mex_bin(...
%     fp_in, fp_out, theta_hat, kappa_hat, gamma, uv_scale_factor);
% 
% muv = read_obj(fp_out);
% V = [muv.vt, zeros(size(muv.vt, 1), 1)]; F = muv.f.vt;
% write_obj(out_folder, [meshname '_ourcycles_ioq_uv' outext], V, F);
% 
% 
% %%
% fp_out = fullfile(out_folder, [meshname '_marcelcycles_ioq' outext]);
% [u2, v2] = conformal_seamless_mex_bin(...
%     fp_in, fp_out, theta_hat, [], [], uv_scale_factor);
% 
% muv = read_obj(fp_out);
% V = [muv.vt, zeros(size(muv.vt, 1), 1)]; F = muv.f.vt;
% write_obj(out_folder, [meshname '_marcelcycles_ioq_uv' outext], V, F);
% 
% %%
% fp_out = fullfile(out_folder, [meshname '_ourcycles_marcelkappa_ioq' outext]);
% [u3, v3] = conformal_seamless_mex_bin(...
%     fp_in, fp_out, theta_hat, [], gamma, uv_scale_factor);
% 
% muv = read_obj(fp_out);
% V = [muv.vt, zeros(size(muv.vt, 1), 1)]; F = muv.f.vt;
% write_obj(out_folder, [meshname '_ourcycles_marcelkappa_ioq_uv' outext], V, F);
% 
% %% MIQ
% theta_hat = (4 - kmiq) * (pi / 2);
% kappa_hat = round(flatten_generators(m) / (pi/2)) * pi/2;
% fp_out = fullfile(out_folder, [meshname '_ourcycles_miq' outext]);
% [u4, v4] = conformal_seamless_mex_bin(...
%     fp_in, fp_out, theta_hat, kappa_hat, gamma, uv_scale_factor);
% 
% muv = read_obj(fp_out);
% V = [muv.vt, zeros(size(muv.vt, 1), 1)]; F = muv.f.vt;
% write_obj(out_folder, [meshname '_ourcycles_miq_uv' outext], V, F);
% 
% %%
% fp_out = fullfile(out_folder, [meshname '_marcelcycles_miq' outext]);
% [u5, v5] = conformal_seamless_mex_bin(...
%     fp_in, fp_out, theta_hat, [], [], uv_scale_factor);
% 
% muv = read_obj(fp_out);
% V = [muv.vt, zeros(size(muv.vt, 1), 1)]; F = muv.f.vt;
% write_obj(out_folder, [meshname '_marcelcycles_miq_uv' outext], V, F);
% 
% %%
% fp_out = fullfile(out_folder, [meshname '_ourcycles_marcelkappa_miq' outext]);
% [u3, v3] = conformal_seamless_mex_bin(...
%     fp_in, fp_out, theta_hat, [], gamma, uv_scale_factor);
% 
% muv = read_obj(fp_out);
% V = [muv.vt, zeros(size(muv.vt, 1), 1)]; F = muv.f.vt;
% write_obj(out_folder, [meshname '_marcelcycles_miq_uv' outext], V, F);
% 
% clear functions; % Free the mex files for writing (so that it's possible to recompile)