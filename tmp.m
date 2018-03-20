% % %%
% % m = Mesh();
% % m.loadTM('../data/bunnies/bunny_99k_faces.off');
% % L = -cotmatrix(m.V, m.F);
% % lb = 8000; % block size
% % A = single(L+1/m.nV);
% % 
% % %%
% % Linv = block_inv_gpu(A, lb);
% % 
% % %%
% % disp('cpu')
% % timeit(@() inv(A))
% % 
% % %% Find best block size
% % for lb = 1000:1000:20000
% %     disp(lb)
% %     T = gputimeit(@() block_inv_gpu(A, lb));
% %     fprintf("blksize: %d, time: %g\n", T);
% % end
% % 
% % %% run OTC
% % fp = '../data/bunnies/bunny_res7.off';
% % g_constraint_vec = [1, 0, 0];
% % opt_block_inv_gpu = struct('block_inv_gpu', true, ...
% %         'CUDAKernel_reduce_cols', true, ...
% %         'TCODS', true, ...
% %         'gConstraintVec', g_constraint_vec, ...
% %         'lb', 20000);
% % f0 = [1]; theta0 = 0; degree = 4; n_iter=1000;
% % m = Mesh();
% % m.loadTM(fp);
% % res_otc = run_lattice_iter(m, ...
% %     f0, ...
% %     theta0, ...
% %     degree, ...
% %     n_iter, ...
% %     opt_block_inv_gpu);
% % res_miq = nrosy_mex(fp, f0, g_constraint_vec, degree);
% % 
% % figure(1)
% % MeshVis.plot(m, ...
% %     'nrosy', res_otc, ...
% %     'ConstrainedFaces', f0, ...
% %     'ConstraintVectors', g_constraint_vec);
% % figure(2)
% % MeshVis.plot(m, ...
% %     'nrosy', res_miq, ...
% %     'ConstrainedFaces', f0, ...
% %     'ConstraintVectors', g_constraint_vec);
% % link_figures(1, 2);
% % 
% % 
% % 
% % 
% % % %% test block inv with block mult
% % % L = -cotmatrix(m.V, m.F);
% % % lb = 20000; % block size
% % % A = single(L+1/m.nV);
% % % Lp1 = block_inv_gpu(A, lb, true);
% % % Lp2 = block_inv_gpu(A, lb, false);
% % % check_norm('Lp1', 'Lp2')
% % 
% % %% Genus tests
% % m = Mesh();
% % m.loadTM('../data/torus_fat_r2.off');
% % assert(abs(m.genus - 1) < 1e-10)
% % m.loadTM('../data/simple/dancer2.off')
% % assert(abs(m.genus - 1) < 1e-10)
% % 
% % %% Filter genus > 0 meshes
% % filepaths = get_filepaths('../data/inputmodels', 'obj');
% % %filepaths = {'..\data\inputmodels\woodenfish.obj'};
% % %filepaths = {'..\data\inputmodels\armadillo.obj'};
% % %filepaths = {'..\data\inputmodels\vh_skin.obj'};
% % is_genus0 = @(m) abs(m.genus) < 1e-10;
% % simple_filepaths = {};
% % for i = 1:numel(filepaths)
% %     fp = filepaths{i};
% %     disp(['Loading ', fp, '...'])
% %     mesh = Mesh();
% %     mesh.loadTM(fp);
% %     try
% %         if is_genus0(mesh)           
% %             [~, name, ext] = fileparts(fp);
% %             out_fp = fullfile('..', 'data', 'genus0', [name, '.off']);
% %             fprintf('\tSaving to %s\n', out_fp);
% %             mesh.saveTM(out_fp);
% %         end
% %     catch ex
% %         fprintf('Caught exception while filtering: %s\r\n', ex.message);
% %     end
% % end
% % 
% % %% Save ffield test
% % m = Mesh();
% % %m.loadTM('../data/bunny.off');
% % m.loadTM('../data/genus0/armadillo.off');
% % S = [1, 0.5; 100, 1.5];
% % f0 = [1];
% % theta0 = 0;
% % degree = 4;
% % res_otc = TCODS(m, S, f0, theta0, degree);
% % 
% % m.set_ffield(...
% %              res_otc.degree, ...
% %              res_otc.local_frames, ...
% %              res_otc.frame_diffs, ...
% %              res_otc.theta, ...
% %              res_otc.ffield, ...
% %              res_otc.S, ...
% %              res_otc.E)
% % 
% % m.saveTM('bunny.ffield');
% % 
% % m2 = Mesh();
% % m2.loadTM('bunny.ffield');
% % 
% % check_norm('m.degree', 'm2.degree', 'Tol', 1e-8);
% % check_norm('m.local_frames', 'm2.local_frames', 'Tol', 1e-8);
% % check_norm('m.frame_diffs', 'm2.frame_diffs', 'Tol', 1e-8);
% % check_norm('m.ffield_angles', 'm2.ffield_angles', 'Tol', 1e-8);
% % check_norm('m.ffield_vectors', 'm2.ffield_vectors', 'Tol', 1e-8);
% % check_norm('m.singularities', 'm2.singularities', 'Tol', 1e-8);
% % check_norm('m.miq_energy', 'm2.miq_energy', 'Tol', 1e-8);
% % 
% % figure()
% % subplot(121)
% % m.draw()
% % title('before save ffield')
% % subplot(122)
% % m2.draw()
% % title('after save ffield')
% % 
% % %%
% % m = Mesh();
% % fp = '../results/experiments/bunnies/ffields/bunny_res7_IOQ_conn.ffield';
% % m.loadTM(fp);
% % V = m.V;
% % F = m.F;
% % R = m.ffield_vectors(1:m.nF, :);
% % [UV, FUV] = MIQ_param_mex_bin(V, F, R);
% % 
% % %%
% % m = Mesh();
% % fp = '../data/bunny.off';
% % m.loadTM(fp);
% % 
% % inv_gpu_range = [-inf, 21582];
% % inv_block_range = [21582, 50002];
% % inv_cpu_range = [50002, inf];
% % iter_gpu_range = [-inf, 21582];
% % alg_selector = create_IOQ_alg_selector(...
% %         inv_gpu_range, ...
% %         inv_block_range, ...
% %         inv_cpu_range, ...
% %         iter_gpu_range);
% % G_CONSTRAINT_VEC = [1, 0, 0];
% % LB = 20000;
% % opt_conn = struct('TCODS', true, ...
% %     'gConstraintVec', G_CONSTRAINT_VEC, ...
% %     'graph', true, ...
% %     'alg_selector', alg_selector, ...
% %     'lb', LB);
% % 
% % opt_cot = struct('TCODS', true, ...
% %     'gConstraintVec', G_CONSTRAINT_VEC, ...
% %     'graph', false, ...
% %     'alg_selector', alg_selector, ...
% %     'lb', LB);
% % 
% % f0 = [1];
% % theta0 = 0;
% % degree = 4;
% % n_iter = 1000;
% % 
% % %res_conn = run_lattice_iter(m, f0, theta0, degree, n_iter, opt_conn);
% % res_cot = run_lattice_iter(m, f0, theta0, degree, n_iter, opt_cot);
% % 
% % %res_conn.E
% % res_cot.E
% % 
% % %% time inv
% % m = Mesh('../data/bunny.off');
% % gd = gpuDevice();
% % 
% % [~, E] = grad(m);
% % Gf = m.Gf;
% % Gv = m.Gv;
% % 
% % tic
% % Lp1 = 0.25 * Gf * (E' \ Gv);
% % toc
% % 
% % L = gpuArray(full(-cotmatrix(m.V, m.F)));
% % tic
% % %Lp2 = 0.25 * Gf * (inv(gpuArray(full(E'))) * Gv);
% % Lp2 = inv(L + 1/m.nV);
% % wait(gd); toc
% % 
% % L = gather(L);
% % tic
% % %Lp3 = 0.25 * inv(Gf * E' * Gv);
% % Lp3 = inv(L + 1/m.nV);
% % toc
% % 
% % tic
% % %Lp4 = 0.25 * Gf * (invChol_mex(E') * Gv);
% % Lp4 = invChol_mex(L + 1/m.nV);
% % toc
% % 
% % %%
% % experiments_setup;
% % meshname = 'armadillo';
% % [V, F] = read_off([data_folder meshname '.off']);
% % nF = size(F, 1);
% % fid = fopen([out_folder_nrosy meshname '_IOQ_conn.rosy']);
% % str = fgetl(fid);
% % [R, Rcount] = fscanf(fid, '%g %g %g', [3 nF]); 
% % R = R';
% % assert(Rcount == 3*nF)
% % fclose(fid);
% % 
% % [uv, fuv] = MIQ_param_mex_bin(V, F, R, 75);
% %     
% % disp('Saving...')
% % %options.nm_file = ['ext/rosy2gridparam' 'cross.png'];
% % options.nm_file = './cross.png';
% % options.face_texcorrd = fuv;
% % options.object_texture = uv;
% % 
% % write_obj_for_quad(out_folder_gridparams, ...
% %     [meshname, '_IOQ_conn', obj_ext], ...
% %     V, ...
% %     F, ...
% %     options);
% % 
% % in_obj = [out_folder_gridparams, meshname, '_IOQ_conn', obj_ext];
% % out_obj = [out_folder_quads, meshname, '_IOQ_conn', obj_ext];
% % system([QEXBIN, ' ', in_obj, ' ', out_obj]);
% % 
% % 
% % %%
% % files = get_filepaths('../data/bunnies');
% % for i = 1:numel(files)
% %     fp = files{i};
% %     disp(['Loading ', fp, '...'])
% %     mesh = Mesh(fp);
% %     L = -cotmatrix(mesh.V, mesh.F);
% %     
% %     disp('pseudoinv')
% %     %timeit(@() pseudoinverse(L))
% %     tic; pseudoinverse(L); toc
% %     
% %     disp('chol')
% %     L = full(L);
% %     %timeit(@() invChol_mex(L+1/nv)-1/nv)
% %     tic; invChol_mex(L+1/nv)-1/nv; toc
% % end
% % 
% % %%
% % base_folder = '../results/experiments/genus0/ffields';
% % meshname = 'dente';
% % fp1 = fullfile(base_folder, [meshname, '_IOQ_conn.ffield']);
% % fp2 = fullfile(base_folder, [meshname, '_MIQ.ffield']);
% % mesh1 = Mesh(fp1);
% % mesh2 = Mesh(fp2);
% % DEGREE = 4;
% % 
% % Kg = get_gaussian_curvature(mesh1);
% % [d0, d1] = get_exterior_derivatives(mesh1);
% % L = d0'*d0;
% % %tic; Lp = pseudoinverse(L); toc
% % 
% % nv1 = mesh1.nV; nf1 = mesh1.nF;
% % sing1 = mesh1.vert_sing;
% % ns1 = size(sing1, 1);
% % k1 = full(sparse(sing1(:, 1), ones(ns1, 1), DEGREE*sing1(:, 2), nv1, 1));
% % b1 = -Kg + (pi / 2) * k1;
% % u1 = L \ b1;
% % u1 = u1 - mean(u1);
% % 
% % check_norm('L*u1', 'b1', 'Log', -1);
% % 
% % nv2 = mesh2.nV; nf2 = mesh2.nF;
% % sing2 = mesh2.vert_sing;
% % ns2 = size(sing2, 1);
% % k2 = full(sparse(sing2(:, 1), ones(ns2, 1), DEGREE*sing2(:, 2), nv2, 1));
% % b2 = -Kg + (pi / 2) * k2;
% % u2 = L \ b2;
% % u2 = u2 - mean(u2);
% % 
% % check_norm('L*u1', 'b1', 'Log', -1);
% % 
% % fprintf('max(u1) - min(u1) = %g\n', max(u1) - min(u1));
% % fprintf('max(u2) - min(u2) = %g\n', max(u2) - min(u2));
% % 
% % %%
% % base_folder = '../results/experiments/genus0/ffields';
% % files = get_filepaths(base_folder);
% % processed = containers.Map();
% % IOQ_maxu_minus_minu = [];
% % MIQ_maxu_minus_minu = [];
% % labels = {};
% % for i = 1:numel(files)
% %     fp = files{i};
% %     disp([fp, '...'])
% %     [~, full_meshname, ~] = fileparts(fp);
% %     idx = strfind(full_meshname, '_');
% %     meshname = full_meshname(1:idx(1)-1);
% %     j = 1;
% %     while ~exist(fullfile(base_folder, [meshname, '_MIQ.ffield']), 'file')
% %         if j > length(idx)
% %             warning(['Could not parse ', fp])
% %             continue
% %         end
% %         meshname = full_meshname(1:idx(j)-1);
% %         j = j + 1;
% %     end
% %     if processed.isKey(meshname) && processed(meshname) == true
% %         continue
% %     end
% %     
% %     fp1 = fullfile(base_folder, [meshname, '_IOQ_conn.ffield']);
% %     fp2 = fullfile(base_folder, [meshname, '_MIQ.ffield']);
% %     try
% %         mesh1 = Mesh(fp1);
% %         mesh2 = Mesh(fp2);
% %     catch
% %         disp('Could not load ')
% %         disp(fp1)
% %         disp(fp2)
% %     end
% %     DEGREE = 4;
% % 
% %     Kg = get_gaussian_curvature(mesh1);
% %     [d0, d1] = get_exterior_derivatives(mesh1);
% %     L = d0'*d0;
% %     %tic; Lp = pseudoinverse(L); toc
% % 
% %     nv1 = mesh1.nV; nf1 = mesh1.nF;
% %     sing1 = mesh1.vert_sing;
% %     ns1 = size(sing1, 1);
% %     k1 = full(sparse(sing1(:, 1), ones(ns1, 1), DEGREE*sing1(:, 2), nv1, 1));
% %     b1 = -Kg + (pi / 2) * k1;
% %     u1 = L \ b1;
% %     u1 = u1 - mean(u1);
% % 
% %     check_norm('L*u1', 'b1', 'Log', -1);
% % 
% %     nv2 = mesh2.nV; nf2 = mesh2.nF;
% %     sing2 = mesh2.vert_sing;
% %     ns2 = size(sing2, 1);
% %     k2 = full(sparse(sing2(:, 1), ones(ns2, 1), DEGREE*sing2(:, 2), nv2, 1));
% %     b2 = -Kg + (pi / 2) * k2;
% %     u2 = L \ b2;
% %     u2 = u2 - mean(u2);
% % 
% %     check_norm('L*u2', 'b2', 'Log', -1);
% % 
% %     fprintf('max(u1) - min(u1) = %g\n', max(u1) - min(u1));
% %     fprintf('max(u2) - min(u2) = %g\n', max(u2) - min(u2));
% %     IOQ_maxu_minus_minu(end+1) = max(u1) - min(u1);
% %     MIQ_maxu_minus_minu(end+1) = max(u2) - min(u2);
% %     
% %     processed(meshname) = true;
% %     labels{end+1} = meshname;
% % end
% % 
% % figure()
% % xx = 1:length(IOQ_maxu_minus_minu);
% % plot(xx, IOQ_maxu_minus_minu, 'xb', 'XTickLabel', labels, ...
% %      xx, MIQ_maxu_minus_minu, 'xr', 'XTickLabel', labels)
% % legend('IOQ max(u) - min(u)', 'MIQ max(u) - min(u)')
% % 
% % figure()
% % plot(xx, IOQ_maxu_minus_minu ./ MIQ_maxu_minus_minu, '-xr')
% % title('[IOQ max(u) - min(u)] / [MIQ max(u) - min(u)]')
% % 
% % %%
% % m = Mesh('../data/torus_s0.off');
% % V = m.V;
% % F = m.F;
% % nv = m.nV; nf = m.nF;
% % cycles = m.generator_cycles();
% % ng2 = numel(cycles);
% % 
% % k = IOC_cpu(V, F);
% % 
% % [A, K, d0, d1, H] = tcods_gsystem(V, F);
% % theta_hat = -K(1:nv) + (pi / 2) * k;
% % kappa_hat = zeros(ng2, 1);
% % 
% % % Create a vector of cycles (face-id based), separated by -1
% % gamma = [];
% % for i = 1:ng2-1
% %     cy = cycles{i};
% %     gamma = [gamma; cy(:); -1];
% % end
% % cy = cycles{end};
% % gamma = [gamma; cy(:)];
% % 
% % conformal_seamless_mex_bin(V, F, theta_hat, kappa_hat, gamma);
% % 
% % %%
% % m = Mesh('../data/torus_s0.off');
% % V = m.V;
% % F = m.F;
% % [A, K, d0, d1, H] = tcods_gsystem(V, F);
% % B = H - d0*(d0\H);
% % Lp = pinv(full(d0'*d0));
% % X = H'*d0*Lp;
% 
% %%
% % FACE0 = 1;        % starting face for tcods
% % THETA0 = nan;     % (use global constraint instead)
% % GVEC = [1, 0, 0]; % constraint vector (in global coordinates)
% % DEGREE = 4;       % works only for degree 4 atm
% % 
% % m = Mesh('../data/elephant_r.off');
% % V = m.V; F = m.F; ng2 = 2*m.genus; nv = m.nV;
% % [A, K, d0, d1, H] = tcods_gsystem(V, F);
% % alpha_G = K(1:end-ng2);
% % beta_G = K(end-ng2+1:end);
% % 
% % rng(112);
% % [alpha_P1,Lp1,E_hist1,m_hist1] = IOQ_cpu_genus0(V, F, 'Laplacian', 'cot','Plot',true,'Iterations',1000);
% % k1 = [alpha_P1; zeros(ng2, 1)];
% % S1v = find(alpha_P1); S1v = [S1v, alpha_P1(S1v)];
% % S1g = [(1:ng2)', zeros(ng2, 1)];
% % 
% % rng(112);
% % [alpha_P2,beta_P2,Lp2,E_hist2,m_hist2] = IOQ_highgenus(V, F, 'Laplacian', 'cot','Plot',true,'Iterations',1000, 'highg_method', 'round');
% % k2 = [alpha_P2; beta_P2];
% % S2v = find(alpha_P2); S2v = [S2v, alpha_P2(S2v)];
% % S2g = find(beta_P2); S2g = [S2g, beta_P2(S2g)];
% % check_norm('alpha_P1','alpha_P2', 'Log', -1);
% % 
% % rng(112);
% % opt.cot = true;
% % res3=run_lattice_iter(m,FACE0,THETA0,DEGREE,1000,opt);
% % k3 = [res3.k; zeros(ng2, 1)];
% % 
% % rng(112);
% % [alpha_P4,Lp4,E_hist4] = IOQ_cpu(V, F, 'Laplacian', 'cot', 'Plot', true);
% % a4 = Lp4 * (alpha_G - (pi/2)*alpha_P4);
% % beta_P4 = round((beta_G - H'*d0*a4) / (pi/2));
% % k4 = [alpha_P4; beta_P4];
% % 
% % %check_norm('Lp1','Lp2');
% % %check_norm('E_hist1(:)','E_hist2(:)');
% % %check_norm('m_hist1(:)','m_hist2(:)');
% % 
% % % [A, K, d0, d1, H] = tcods_gsystem(V, F);
% % % 
% % % n = size(A, 2); 
% % % b1 = K - (pi / 2) * k1;
% % % x1 = lsqlin(speye(n, n), zeros(n, 1), [], [], A, -b1, -inf(n, 1), inf(n, 1));
% % % E1 = norm(x1)^2;
% % % assert(norm(A*x1 - (-b1)) < 1e-10);
% % % assert(norm(d1*x1) < 1e-10);
% % % 
% % % b2 = K - (pi / 2) * k2;
% % % x2 = lsqlin(speye(n, n), zeros(n, 1), [], [], A, -b2, -inf(n, 1), inf(n, 1));
% % % E2 = norm(x2)^2;
% % % assert(norm(A*x2 - (-b2)) < 1e-10);
% % % assert(norm(d1*x2) < 1e-10);
% % 
% % %[m1] = TCODS(m, S1v, S1g, FACE0, THETA0, DEGREE, 'gConstraintVec', GVEC);
% % m1 = TCODS(m, 'k', k1, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE);
% % m2 = TCODS(m, 'k', k2, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE);
% % m3 = TCODS(m, 'k', k3, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE);
% % m4 = TCODS(m, 'k', k4, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE);
% % fprintf('E1 : %g\nE2 : %g\nE3 : %g\nE4 : %g\n\n', m1.miq_energy, m2.miq_energy, m3.miq_energy, m4.miq_energy);
% % figure()
% % subplot(221); m1.draw; title(['E1 = ', num2str(m1.miq_energy)]);
% % subplot(222); m2.draw; title(['E2 = ', num2str(m2.miq_energy)]);
% % subplot(223); m3.draw; title(['E3 = ', num2str(m3.miq_energy)]);
% % subplot(224); m4.draw; title(['E4 = ', num2str(m4.miq_energy)]);
% % 
% % %%
% % n = 4;
% % A = magic(n);
% % COL_TITLES = {'col1', 'col2', 'col3', 'col4'};
% % ROW_TITLES = {'row1', 'row2', 'row3', 'row4'};
% % T = cell(n+1, n+1);
% % T(1, 2:end) = COL_TITLES;
% % T(2:end, 1) = ROW_TITLES;
% % T(2:end, 2:end) = arrayfun(@(x) num2str(x), A, 'UniformOutput', false);
% % [~, bold] = min(A, [], 2);
% % for i = 1:n
% %     j = bold(i);
% %     T{i+1, j+1} = ['$\textbf{', T{i+1,j+1}, '}$'];
% % end
% % 
% % clear input;
% % input.data = table(T);
% % input.tableCaption = caption_T;
% % input.tablePlacement = 'H';
% % input.tableColumnAlignment = 'c';
% % latex = latexTable(input);
% % 
% 
% %%
% BLOCK_SIZE = 8000;
% LAP_TYPE = 'conn';
% PLOT = false;
% FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0];
% INIT_BETA_P = 'round';
% N_ITER = 1000;
% SEED = 112;
% fp = '../data/genus1_small/3holes.off';
% m = Mesh(fp); nv = m.nV; V = m.V; F = m.F; ng2 = m.genus*2;
% [d0, ~] = get_exterior_derivatives(m);
% L = d0'*d0;
% 
% tic
% Lp = double(block_inv_gpu(full(L + 1/nv), BLOCK_SIZE) - 1/nv);
% toc
% 
% print_header('1');
% rng(SEED);
% [alpha_P1,~,E_hist1,m_hist1] = IOQ_genus0_gpu(V, F, 'LaplacianPinv', Lp,'Plot',PLOT,'Iterations',N_ITER);
% k1 = [alpha_P1; zeros(ng2, 1)];
% m1 = TCODS(m, 'k', k1, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', true);
% Eioq1 = m1.miq_energy
% 
% print_header('2');
% rng(SEED);
% [alpha_P2,beta_P2,~,E_hist2,m_hist2] = IOQ_highgenus_gpu(V, F, 'LaplacianPinv', Lp,'Plot',PLOT,'Iterations',N_ITER, 'highg_method', 'option1a', 'beta_P', INIT_BETA_P);
% k2 = [alpha_P2; beta_P2];
% m2 = TCODS(m, 'k', k2, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', true);
% Eioq2 = m2.miq_energy;
% 
% print_header('3');
% rng(SEED);
% [alpha_P3,beta_P3,~,E_hist3,m_hist3] = IOQ_highgenus_gpu(V, F, 'LaplacianPinv', Lp,'Plot',PLOT,'Iterations',N_ITER, 'highg_method', 'option1b', 'beta_P', INIT_BETA_P);
% k3 = [alpha_P3; beta_P3];
% m3 = TCODS(m, 'k', k3, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', true);
% Eioq3 = m3.miq_energy;
% 
% print_header('4');
% rng(SEED);
% [alpha_P4,beta_P4,~,E_hist4,m_hist4] = IOQ_highgenus_gpu(V, F, 'LaplacianPinv', Lp,'Plot',PLOT,'Iterations',N_ITER, 'highg_method', 'option2_optimized', 'beta_P', INIT_BETA_P);
% k4 = [alpha_P4; beta_P4];
% m4 = TCODS(m, 'k', k4, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', true);
% Eioq4 = m4.miq_energy;
% 
% print_header('miq');
% [theta, p, kmiq, R, local_frames, frame_diffs, elapsed, pFixed] = ...
%         NRoSy_mex_bin(fp, FACE0, GVEC, DEGREE);
% Emiq = E_MIQ(m, theta, frame_diffs, p, DEGREE);
% 
% E = [Eioq1; Eioq2; Eioq3; Eioq4; Emiq]
% 
% %%
% % filepaths = get_filepaths('../data/genus1_small', '.off');
% % for i = 1:numel(filepaths)
% %     fp = filepaths{i};
% %     print_header(fp);
% %     mesh = Mesh(fp);
% %     flatten_generators(mesh);
% % end
% 
% fp = '../data/genus1_small/bucky180.off';
% m = Mesh(fp);
% flatten_generators(m);
% cycles = m.generator_cycles;
% ng2 = 2*m.genus;
% FE = m.FEAdj;
% %%
% for i = 1:ng2
%     disp(i)
%     cy = cycles{i};
%     for j = 1:length(cy)-1
%         f1 = cy(j);
%         f2 = cy(j+1);
%         e1 = FE(f1, :);
%         e2 = FE(f2, :);
%         common_edge = intersect(e1, e2);
%         assert(length(common_edge)==1)
%     end
% end
% 
% figure(); m.colorFaces(cycles{1}, 'r');
% 
% %% torus
% BLOCK_SIZE = 8000;
% LAP_TYPE = 'conn';
% PLOT = true;
% FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0];
% INIT_BETA_P = 'round';
% N_ITER = 1000;
% SEED = 112;
% fp = '../data/genus1_small/torus_s0.off';
% m = Mesh(fp); nv = m.nV; V = m.V; F = m.F; ng2 = m.genus*2;
% [d0, ~] = get_exterior_derivatives(m);
% L = d0'*d0;
% 
% gd = gpuDevice();
% tic
% Lp = double(block_inv_gpu(full(L + 1/nv), BLOCK_SIZE) - 1/nv);
% wait(gd); toc
% %Lp = inv(L+1/nv) - 1/nv;
% 
% print_header('1');
% rng(SEED);
% [alpha_P1,beta_P,~,E_hist1,m_hist1] = IOQ_highgenus(V, F, 'LaplacianPinv', Lp,'Plot',PLOT,'Iterations',N_ITER, 'highg_method', 'option1a');
% k1 = [alpha_P1; zeros(ng2, 1)];
% m1 = TCODS(m, 'k', k1, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', true);
% Eioq1 = m1.miq_energy
% 
% %% 3holes
% BLOCK_SIZE = 8000;
% LAP_TYPE = 'conn';
% PLOT = false;
% FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0];
% INIT_BETA_P = 'round';
% %INIT_BETA_P = [];
% N_ITER = 1000;
% SEED = 112;
% fp = '../data/bucky180_r.off';
% m = Mesh(fp); nv = m.nV; V = m.V; F = m.F; ng2 = m.genus*2;
% %[d0, ~] = get_exterior_derivatives(m);
% %L = d0'*d0;
% 
% %[alpha_P2,beta_P2,Lp2,E_hist2,m_hist2] = ...
% %    IOQ_highgenus(V, F, 'Laplacian', LAP_TYPE,'Plot',PLOT,'Iterations',1000, 'highg_method', 'option1a', 'beta_P', INIT_BETA_P);
% %tic
% %Lp = double(block_inv_gpu(full(L + 1/nv), BLOCK_SIZE) - 1/nv);
% %toc
% %Lp = inv(L+1/nv) - 1/nv;
% 
% print_header('1');
% rng(SEED);
% [alpha_P1,beta_P,~,E_hist1,m_hist1] = IOQ_highgenus_gpu(...
%     V, F, ...
%     'Laplacian', LAP_TYPE, ...
%     'Plot', PLOT, ...
%     'Iterations', N_ITER, ...
%     'highg_method', 'option1a', ...
%     'beta_P', INIT_BETA_P, ...
%     'n_alternating', 10);
% k1 = [alpha_P1; beta_P];
% m1 = TCODS(m, 'k', k1, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
% Eioq1 = m1.miq_energy
% 
% print_header('miq');
% [theta, p, kmiq, R, local_frames, frame_diffs, elapsed, pFixed] = ...
%         NRoSy_mex_bin(fp, FACE0, GVEC, DEGREE);
% Emiq = E_MIQ(m, theta, frame_diffs, p, DEGREE);
% 
% fprintf('Eioq1 = %g\nEmiq = %g\n', Eioq1, Emiq);
% 
% %% Compare cpu to gpu
% BLOCK_SIZE = 8000;
% LAP_TYPE = 'conn';
% PLOT = true;
% FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0];
% %INIT_BETA_P = 'round';
% INIT_BETA_P = [];
% N_ITER = 1000;
% SEED = 112;
% fp = '../data/genus1_small/3holes.off';
% m = Mesh(fp); nv = m.nV; V = m.V; F = m.F; ng2 = m.genus*2;
% [d0, ~] = get_exterior_derivatives(m);
% 
% print_header('cpu')
% [alpha_P1,beta_P1,~,E_hist1,m_hist1] = IOQ_highgenus(...
%     V, F, ...
%     'Laplacian', LAP_TYPE, ...
%     'Plot', PLOT, ...
%     'Iterations', N_ITER, ...
%     'highg_method', 'option1a', ...
%     'beta_P', INIT_BETA_P);
% k1 = [alpha_P1; beta_P1];
% m1 = TCODS(m, 'k', k1, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
% Eioq1 = m1.miq_energy;
% 
% 
% print_header('gpu')
% [alpha_P2, beta_P2, ~, E_hist2, m_hist2] = IOQ_highgenus_gpu(...
%     V, F, ...
%     'Laplacian', LAP_TYPE, ...
%     'Plot', PLOT, ...
%     'Iterations', N_ITER, ...
%     'highg_method', 'option1a', ...
%     'beta_P', INIT_BETA_P);
% k2 = [alpha_P2; beta_P2];
% m2 = TCODS(m, 'k', k2, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
% Eioq2 = m2.miq_energy;
% 
% fprintf('E cpu = %g\n', Eioq1);
% fprintf('E cpu = %g\n', Eioq2);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% test if the results are the same with cot/conn on remeshesd model
% % 5.12.17: seems that no...
% fp = '../data/genus1_small/3holes_r.off';
% FACE0=1;THETA0=0;DEGREE=4;GVEC=[1,0,0];
% m = Mesh(fp); V = m.V; F = m.F;
% method = 'option1a';
% rng(SEED); 
% [alpha_p1, beta_p] = IOQ_highgenus_gpu(...
%     V, F, ...
%     'Laplacian', 'conn', ...
%     'Plot', false, ...
%     'highg_method', method, ...
%     'beta_P', 'round');
% k1 = [alpha_p1; beta_p];
% res1 = TCODS(m, 'k', k1, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
% Eioq1 = res1.miq_energy;
% 
% rng(SEED); 
% [alpha_p2, beta_p2] = IOQ_highgenus_gpu(...
%     V, F, ...
%     'Laplacian', 'cot', ...
%     'Plot', false, ...
%     'highg_method', method, ...
%     'beta_P', 'round');
% k2 = [alpha_p2; beta_p2];
% res2 = TCODS(m, 'k', k2, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
% Eioq2 = res2.miq_energy;
% 
% [theta, p, kmiq, R, local_frames, frame_diffs, elapsed, pFixed] = ...
%         NRoSy_mex_bin(fp, FACE0, GVEC, DEGREE);
% Emiq = E_MIQ(m, theta, frame_diffs, p, DEGREE);
% 
% %fprintf('Eioq1 = %g\nEmiq = %g\n', Eioq1, Emiq);
% 
% %title('Option 2 (cpu, conn), 3holes_r')
% 
% fprintf('Eioq1 = %g\nEioq2 = %g\nEmiq = %g\n', Eioq1, Eioq2, Emiq);
% 
% %title('Option 2 (cpu, cot), 3holes_r')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% opt 1 vs opt 2
% % trying to understand why optimizing E1 yields better results than
% % optimizing E1+E2
% %fp = '../data/genus1_small/3holes.off';
% %fp = '../data/genus1_small/torus_s0.off';
% fp = '../data/bucky_r_1000v.off';
% m = Mesh(fp); V = m.V; F = m.F;
% method = 'option2';
% SEED = 112;FACE0=1;THETA0=0;DEGREE=4;GVEC=[1,0,0];
% rng(SEED); 
% [alpha_p1, beta_p1, ~, out1] = IOQ_highgenus_gpu(...
%     V, F, ...
%     'Laplacian', 'conn', ...
%     'Plot', true, ...
%     'highg_method', 'option1a', ...
%     'beta_P', 'round', ...
%     'Debug', false, ...
%     'UseGPU', false, ...
%     'InvMethod', 'CholMexInv');
% k1 = [alpha_p1; beta_p1];
% res1 = TCODS(m, 'k', k1, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
% Eioq1 = res1.miq_energy;
% 
% title('Option 2 (cpu, conn), bucky_r100v')
% 
% rng(SEED); 
% [alpha_p2, beta_p2, ~, out2] = IOQ_highgenus_gpu(...
%     V, F, ...
%     'Laplacian', 'conn', ...
%     'Plot', true, ...
%     'highg_method', 'option2', ...
%     'beta_P', 'round', ...
%     'n_alternating', 10, ...
%     'Debug', false, ...
%     'UseGPU', false, ...
%     'InvMethod', 'CholMexInv');
% k2 = [alpha_p2; beta_p2];
% res2 = TCODS(m, 'k', k2, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
% Eioq2 = res2.miq_energy;
% 
% figure
% subplot(121)
% E_hist1 = out1.E_hist1;
% E_hist2 = out1.E_hist2;
% E_round_hist2 = out1.E_round_hist2;
% Emiq_hist = out1.Emiq_hist;
% Emiq_round_hist = out1.Emiq_round_hist;
% m_hist = out1.m_hist;
% gridsearch_ticks = out1.gridsearch_ticks;
% xx = 1:length(E_hist1);
% hold on
% lgd = {};
% if ~isempty(E_hist1) && ~isempty(E_hist2)
%     plot(xx, E_hist1+E_hist2, 'LineWidth', 2)
%     lgd{end+1} = 'E1+E2';
% end
% if ~isempty(E_hist1)
%     plot(xx, E_hist1, 'LineWidth', 2)
%     lgd{end+1} = 'E1';
% end
% if ~isempty(E_hist2)
%     plot(xx, E_hist2, 'LineWidth', 2)
%     lgd{end+1} = 'E2';
% end
% if ~isempty(E_round_hist2)
%     plot(xx, E_round_hist2, 'LineWidth', 2)
%     lgd{end+1} = 'E2 w/ y=round(...)';
% end
% if ~isempty(Emiq_hist)
%     plot(xx, Emiq_hist, 'LineWidth', 2)
%     lgd{end+1} = 'Emiq';
% end
% if ~isempty(Emiq_round_hist)
%     plot(xx, Emiq_round_hist, 'LineWidth', 2)
%     lgd{end+1} = 'Emiq w/ y=round(..)';
% end
% if ~isempty(m_hist)
%     plot(xx, m_hist)
%     lgd{end+1} = 'm';
% end
% yl = ylim();
% xx = find(gridsearch_ticks);
% plot(xx, ones(1, length(xx))*yl(1), 'rX')
% hold off
% legend(lgd{:})
% xlabel('Iteration')
% title(['opt 1, final Emiq = ', num2str(Eioq1)])
% 
% subplot(122)
% %figure
% E_hist1 = out2.E_hist1;
% E_hist2 = out2.E_hist2;
% E_round_hist2 = out2.E_round_hist2;
% Emiq_hist = out2.Emiq_hist;
% Emiq_round_hist = out2.Emiq_round_hist;
% m_hist = out2.m_hist;
% gridsearch_ticks = out2.gridsearch_ticks;
% xx = 1:length(E_hist1);
% hold on
% lgd = {};
% if ~isempty(E_hist1) && ~isempty(E_hist2)
%     plot(xx, E_hist1+E_hist2, 'LineWidth', 2)
%     lgd{end+1} = 'E1+E2';
% end
% if ~isempty(E_hist1)
%     plot(xx, E_hist1, 'LineWidth', 2)
%     lgd{end+1} = 'E1';
% end
% if ~isempty(E_hist2)
%     plot(xx, E_hist2, 'LineWidth', 2)
%     lgd{end+1} = 'E2';
% end
% if ~isempty(E_round_hist2)
%     plot(xx, E_round_hist2, 'LineWidth', 2)
%     lgd{end+1} = 'E2 w/ y=round(...)';
% end
% if ~isempty(Emiq_hist)
%     plot(xx, Emiq_hist, 'LineWidth', 2)
%     lgd{end+1} = 'Emiq';
% end
% if ~isempty(Emiq_round_hist)
%     plot(xx, Emiq_round_hist, 'LineWidth', 2)
%     lgd{end+1} = 'Emiq w/ y=round(..)';
% end
% if ~isempty(m_hist)
%     plot(xx, m_hist)
%     lgd{end+1} = 'm';
% end
% ylim(yl);
% %yl = ylim();
% xx = find(gridsearch_ticks);
% plot(xx, ones(1, length(xx))*yl(1), 'rX')
% hold off
% legend(lgd{:})
% xlabel('Iteration')
% title(['opt 2, final Emiq = ', num2str(Eioq2)])
% 
% 
% %% up to first convergence
% figure
% subplot(121)
% gridsearch_ticks = out1.gridsearch_ticks;
% idx = find(gridsearch_ticks);
% if length(idx) > 1
%     idx = idx(1) : idx(2);
% else
%     idx = 1:length(gridsearch_ticks);
% end
% gridsearch_ticks = gridsearch_ticks(idx);
% E_hist1 = out1.E_hist1(idx);
% E_hist2 = out1.E_hist2(idx);
% E_round_hist2 = out1.E_round_hist2(idx);
% Emiq_hist = out1.Emiq_hist(idx);
% Emiq_round_hist = out1.Emiq_round_hist(idx);
% m_hist = out1.m_hist(idx);
% 
% xx = idx;
% hold on
% lgd = {};
% if ~isempty(E_hist1) && ~isempty(E_hist2)
%     plot(xx, E_hist1+E_hist2, 'LineWidth', 2)
%     lgd{end+1} = 'E1+E2';
% end
% if ~isempty(E_hist1)
%     plot(xx, E_hist1, 'LineWidth', 2)
%     lgd{end+1} = 'E1';
% end
% if ~isempty(E_hist2)
%     plot(xx, E_hist2, 'LineWidth', 2)
%     lgd{end+1} = 'E2';
% end
% if ~isempty(E_round_hist2)
%     plot(xx, E_round_hist2, 'LineWidth', 2)
%     lgd{end+1} = 'E2 w/ y=round(...)';
% end
% if ~isempty(Emiq_hist)
%     plot(xx, Emiq_hist, 'LineWidth', 2)
%     lgd{end+1} = 'Emiq';
% end
% if ~isempty(Emiq_round_hist)
%     plot(xx, Emiq_round_hist, 'LineWidth', 2)
%     lgd{end+1} = 'Emiq w/ y=round(..)';
% end
% if ~isempty(m_hist)
%     plot(xx, m_hist)
%     lgd{end+1} = 'm';
% end
% %yl = ylim();
% %xx = find(gridsearch_ticks);
% %plot(xx, ones(1, length(xx))*yl(1), 'rX')
% hold off
% legend(lgd{:})
% xlabel('Iteration')
% title(['opt 1, final Emiq = ', num2str(Eioq1)])
% %yticks(min(yticks):100:max(yticks));
% %xticks(min(xticks):20:max(xticks));
% 
% subplot(122)
% %figure
% gridsearch_ticks = out2.gridsearch_ticks;
% idx = find(gridsearch_ticks);
% if length(idx) > 1
%     idx = idx(1) : idx(2);
% else
%     idx = 1:length(gridsearch_ticks);
% end
% gridsearch_ticks = gridsearch_ticks(idx);
% E_hist1 = out2.E_hist1(idx);
% E_hist2 = out2.E_hist2(idx);
% E_round_hist2 = out2.E_round_hist2(idx);
% Emiq_hist = out2.Emiq_hist(idx);
% Emiq_round_hist = out2.Emiq_round_hist(idx);
% m_hist = out2.m_hist(idx);
% xx = 1:length(idx);
% hold on
% lgd = {};
% if ~isempty(E_hist1) && ~isempty(E_hist2)
%     plot(xx, E_hist1+E_hist2, 'LineWidth', 2)
%     lgd{end+1} = 'E1+E2';
% end
% if ~isempty(E_hist1)
%     plot(xx, E_hist1, 'LineWidth', 2)
%     lgd{end+1} = 'E1';
% end
% if ~isempty(E_hist2)
%     plot(xx, E_hist2, 'LineWidth', 2)
%     lgd{end+1} = 'E2';
% end
% if ~isempty(E_round_hist2)
%     plot(xx, E_round_hist2, 'LineWidth', 2)
%     lgd{end+1} = 'E2 w/ y=round(...)';
% end
% if ~isempty(Emiq_hist)
%     plot(xx, Emiq_hist, 'LineWidth', 2)
%     lgd{end+1} = 'Emiq';
% end
% if ~isempty(Emiq_round_hist)
%     plot(xx, Emiq_round_hist, 'LineWidth', 2)
%     lgd{end+1} = 'Emiq w/ y=round(..)';
% end
% if ~isempty(m_hist)
%     plot(xx, m_hist)
%     lgd{end+1} = 'm';
% end
% %ylim(yl);
% %yl = ylim();
% %xx = find(gridsearch_ticks);
% %plot(xx, ones(1, length(xx))*yl(1), 'rX')
% hold off
% legend(lgd{:})
% xlabel('Iteration')
% title(['opt 2, final Emiq = ', num2str(Eioq2)])
% %yticks(min(yticks):100:max(yticks));
% %xticks(min(xticks):20:max(xticks));
% 
% %%
% fp = '../data/ashish_nob/blade.off';
% m = Mesh(fp);
% 
% %
% disp('d0, L')
% tic
% d0 = get_exterior_derivatives(m);
% L = d0' * d0;
% toc
% 
% disp('chol')
% tic
% chol(L);
% toc
% 
% disp('gpu chol')
% tic
% chol(gpuArray(single(full(L))));
% toc
% 
% disp('ichol')
% tic
% ichol(L,struct('type','ict','droptol',1e-04,'michol','off'));
% toc
% 
% %% ========================================================================
% %  IOQ with block-wise bsx 
% %% ========================================================================
% FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
% fp = '../data/bunnies/bunny_99k_faces.off';
% m = Mesh(fp);
% V = m.V; F = m.F; nv = m.nV; ne = m.nE;
% d0 = get_exterior_derivatives(m);
% L = d0' * d0;
% %Lp = invChol_mex(L + 1/nv) - 1/nv;
% 
% %
% rng(SEED);
% [alpha_p, beta_p, elapsed_total1] = IOQ_highgenus_gpu(V, F, ...
%                         'highg_method', 'genus0', ...
%                         'Mesh', m, ...
%                         'UseGPU', true, ...
%                         'InvMethod', 'ApproxResistance', ...
%                         'JLEps', 0.5);
% k = [alpha_p; beta_p];
% res1 = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
% 
% rng(SEED);
% [alpha_p, beta_p, elapsed_total2] = IOQ_highgenus_gpu(V, F, ...
%                         'highg_method', 'genus0', ...
%                         'bsx', false, ...
%                         'Mesh', m, ...
%                         'UseGPU', true, ...
%                         'InvMethod', 'CholMexInv', ...
%                         'LaplacianPInv', []);
% k2 = [alpha_p; beta_p];
% res2 = TCODS(m, 'k', k2, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', false);
% [theta, p, kmiq, R, local_frames, frame_diffs, elapsed, pFixed] = ...
%     NRoSy_mex_bin(fp, FACE0, GVEC, DEGREE);
% Emiq = E_MIQ(m, theta, frame_diffs, p, DEGREE);
% 
% figure();
% subplot(121); res1.draw; title(['res1, E = ', num2str(res1.miq_energy)])
% subplot(122); res2.draw; title(['res2, E = ', num2str(res2.miq_energy)])
% 
% elapsed_total1
% elapsed_total2
% res1.miq_energy
% res2.miq_energy
% Emiq
% 
% %% ========================================================================
% %  Test factorize package
% %% ========================================================================
% fp = '../data/bunnies/bunny_99k_faces.off';
% m = Mesh(fp); nv = m.nV; ne = m.nE;
% d0 = get_exterior_derivatives(m);
% L = d0'*d0;
% b = rand(nv, 1);
% b = b - mean(b);
% %
% tic; F = factorize(L); toc
% tic; S = inverse(L)  ; toc
% %
% tic ; F = factorize(L) ; F \ b ; toc
% tic ; S = inverse(L)   ; S * b ; toc
% %%
% timeit(@() F \ b)
% timeit(@() S * b)

%% ========================================================================
%
%  ========================================================================
FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
%fp = '../data/torus_s0.off';
fp = '../data/yeah_right.off';
METHOD = 'option1a';
USE_GPU = true;
%fp = '../data/bunny.off';
disp('Loading...')
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE;
disp('Laplacian...')
d0 = get_exterior_derivatives(m);
L = d0' * d0;
%
%disp('Inverting...')
%tic
%%Lp = inv(L+(1/nv)*speye(nv)) - (1/nv)*speye(nv);
%%Lp = inv(L+1/nv) - 1/nv;
%Lp = inv(gpuArray(L+1/nv))-1/nv;
%elapsed_inv = toc;
load('D:\nahum\laplacian_pinvs\yeah_right_Lp_conn.mat')
Lp = full(double(Lp));

% IOQ 1a
disp('IOQ 1a...')
rng(SEED);
[alpha_p, beta_p, elapsed_total2] = IOQ_highgenus_gpu(V, F, ...
                        'highg_method', METHOD, ...
                        'bsx', false, ...
                        'Mesh', m, ...
                        'UseGPU', USE_GPU, ...
                        'LaplacianPInv', Lp, ...
                        'InvMethod', 'GPUInv', ...
                        'Iterations', 2000);
k = [alpha_p; beta_p];
res2 = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', true, 'Duplicate', true);
[E_edges2, Emiq2] = per_edge_energy(res2);



%% ========================================================================
%  JL IOQ high genus
%  ========================================================================
FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
fp = '../data/torus_s0.off';
%fp = '../data/yeah_right.off';
METHOD = 'option1a';
USE_GPU = false;
%fp = '../data/bunny.off';
disp('Loading...')
m = Mesh(fp);
V = m.V; F = m.F; nv = m.nV; ne = m.nE;
disp('Laplacian...')
d0 = get_exterior_derivatives(m);
L = d0' * d0;
%
disp('Inverting...')
tic
%Lp = inv(L+(1/nv)*speye(nv)) - (1/nv)*speye(nv);
%Lp = inv(L+1/nv) - 1/nv;
Lp = invChol_mex(L+1/nv)-1/nv;
elapsed_inv = toc;

%
% JL IOQ
disp('JL IOQ...')
rng(SEED);
[alpha_p, beta_p, elapsed_total1] = IOQ_highgenus_gpu(V, F, ...
                        'highg_method', METHOD, ...
                        'Mesh', m, ...
                        'UseGPU', USE_GPU, ...
                        'InvMethod', 'ApproxResistance', ...
                        'JLEps', 0.5, ...
                        'Iterations', 2000);
k = [alpha_p; beta_p];
res1 = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', true, 'Duplicate', true);
[E_edges1, Emiq1] = per_edge_energy(res1);

% IOQ 1a
disp('IOQ 1a...')
rng(SEED);
[alpha_p, beta_p, elapsed_total2] = IOQ_highgenus_gpu(V, F, ...
                        'highg_method', 'option1a', ...
                        'bsx', false, ...
                        'Mesh', m, ...
                        'UseGPU', USE_GPU, ...
                        'LaplacianPInv', Lp, ...
                        'InvMethod', 'GPUInv');
k = [alpha_p; beta_p];
res2 = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', true, 'Duplicate', true);
[E_edges2, Emiq2] = per_edge_energy(res2);
elapsed_total2 = elapsed_total2 + elapsed_inv;

% IOQ 2
disp('IOQ 2...')
rng(SEED);
[alpha_p, beta_p, elapsed_total3] = IOQ_highgenus_gpu(V, F, ...
                        'highg_method', 'option2', ...
                        'bsx', false, ...
                        'Mesh', m, ...
                        'UseGPU', true, ...
                        'LaplacianPInv', Lp, ...
                        'InvMethod', 'GPUInv');
k = [alpha_p; beta_p];
res3 = TCODS(m, 'k', k, 'f0', FACE0, 'theta0', THETA0, 'degree', DEGREE, 'CreateFField', true, 'Duplicate', true);
[E_edges3, Emiq3] = per_edge_energy(res3);

% MIQ
disp('MIQ...')
[res4, elapsed_total4] = ...
    nrosy_mex(fp, FACE0, GVEC, DEGREE);
%res3.ffield_vectors = [];
%Emiq = E_MIQ(m, theta, frame_diffs, p, DEGREE);

disp('GO...')
tic
[Ego4, Emiq4, x4, ffield4, S4] = calc_GO_and_MIQ_energies(fp, 4);
elapsed_total5 = toc;
res5 = Mesh(fp);
res5.set_ffield(DEGREE, ...
                [], ...
                [], ...
                [], ...
                ffield4, ...
                S4, ...
                [], ...
                Emiq4, ...
                []);

% Plot
disp('Plotting')
plotting_defaults
title1 = {'JL IOQ', ...
    sprintf('Etc = %.4g', res1.miq_energy), ...
    sprintf('Emiq = %.4g', Emiq1), ...
    sprintf('T = %.4g', elapsed_total1), ...
    sprintf('# sing = %d', res1.n_vert_sing)};
title2 = {'IOQ 1a', ...
    sprintf('Etc = %.4g', res2.miq_energy), ...
    sprintf('Emiq = %.4g', Emiq2), ...
    sprintf('T = %.4g', elapsed_total2), ...
    sprintf('# sing = %d', res2.n_vert_sing)};
title3 = {'IOQ 2', ...
    sprintf('Etc = %.4g', res3.miq_energy), ...
    sprintf('Emiq = %.4g', Emiq3), ...
    sprintf('T = %.4g', elapsed_total3), ...
    sprintf('# sing = %d', res3.n_vert_sing)};
title4 = {'MIQ', ...
    sprintf('Emiq = %.4g', res4.miq_energy), ...
    sprintf('T = %.4g', elapsed_total4), ...
    sprintf('# sing = %d', res4.n_vert_sing)};
title5 = {'GO', ...
    sprintf('Emiq = %.4g', Emiq4), ...
    sprintf('Ego = %.4g', Ego4), ...
    sprintf('T = %.4g', elapsed_total5), ...
    sprintf('# sing = %d', size(S4, 1))};

figure
subplot(242); res1.draw('FaceAlpha', 1); subplot(241); title(title1); axis off
subplot(244); res2.draw('FaceAlpha', 1); subplot(243); title(title2); axis off
subplot(246); res5.draw('FaceAlpha', 1); subplot(245); title(title5); axis off
subplot(248); res4.draw('FaceAlpha', 1); subplot(247); title(title4); axis off

