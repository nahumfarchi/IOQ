function [alg_selector, ...
	inv_gpu_range, inv_block_range, inv_cpu_range, iter_gpu_range] = find_IOQ_ranges(filepaths, lb)

	inv_gpu_range = [-inf, 0];
	inv_block_range = [-inf, nan];
	inv_cpu_range = [inf, inf];
	iter_gpu_range = [-inf, 0];

	n_files = numel(filepaths);
	for i = 1:n_files
		fp = filepaths{i};
        disp(fp)
		m = Mesh();
		m.loadTM(fp);

        try
            L = single(full(-cotmatrix(m.V, m.F)));
        catch
            disp('Too big (skipping)...')
            continue
        end

		try % Invert on gpu
			if ~isnan(inv_block_range(2)) && m.nV > inv_block_range(2)
				error('go to block')
			end
			if m.nV < inv_gpu_range(2)
				disp('skipping...')
				continue
			end

			Lp = inv(gpuArray(L) + 1/m.nV);
            disp('gpu')
			clear Lp;

			inv_gpu_range(2) = max(inv_gpu_range(2), m.nV);
			iter_gpu_range(2) = max(iter_gpu_range(2), m.nV);
		catch
			try % Block inversion
				if m.nV > inv_cpu_range(1)
					error('go to cpu')
				end
				if ~isnan(inv_block_range(2)) && m.nV < inv_block_range(2)
					disp('skipping...')
					continue
				end

				Lp = block_inv_gpu(L + 1/m.nV, lb); % Returns a non-gpu array
                disp('block')
				inv_block_range(2) = max(inv_block_range(2), m.nV);

				try % See if Lp fits on the gpu
					Lp = gpuArray(Lp);
					clear Lp;
					iter_gpu_range(2) = max(iter_gpu_range(2), m.nV);
				catch

				end
			catch % Cpu inversion
                disp('cpu')
				inv_cpu_range(1) = min(inv_cpu_range(1), m.nV);
			end
		end
	end

	inv_block_range(1) = inv_gpu_range(2);
	inv_cpu_range(1) = inv_block_range(2);

    alg_selector = create_IOQ_alg_selector(...
        inv_gpu_range, ...
        inv_block_range, ...
        inv_cpu_range, ...
        iter_gpu_range);

	% jump = 1000;

	% loop_gpu_range = [0, 0];
	% ME = [];
	% sz = 1000;
	% while isempty(ME)
	% 	try
	% 		A = gpuArray(zeros(sz, sz, CLASSNAME));
	% 		clear A;
	% 		sz = sz + jump;
	% 	catch ME
	% 		loop_gpu_range(2) = sz - jump;
	% 	end
	% end

	% inv_gpu_range = [0, 0];
	% ME = [];
	% sz = 1000;
	% while isempty(ME)
	% 	try
	% 		A = gpuArray(zeros(sz, sz, CLASSNAME))
	% 		A = inv(A);
	% 		clear A;
	% 	catch ME
	% 		inv_gpu_range(2) = sz - jump;
	% 	end
	% end

	% inv_block_range = [0, 0];
	% ME = [];
	% sz = 1000;
	% while isempty(ME)
	% 	try
	% 		A = gpuArray(zeros(sz, sz, CLASSNAME));
	% 		A = block_inv_gpu(A)
end