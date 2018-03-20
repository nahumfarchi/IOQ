function [defects, flat_meshes] = flatten_cycles(mesh, cycles, debug)
    %FLATTEN_GENERATORS
    
    if nargin < 3
        debug = false;
    end

    flat_meshes = {};
    ncycles = numel(cycles);
    defects = zeros(ncycles, 1);
    
    for i = 1:ncycles
        V = mesh.V; F = mesh.F;
        nv = mesh.nV; nf = mesh.nF;
        
        cy = cycles{i};
        cy = cy(1:end-1);
        cy_len = length(cy);
        if cy_len <= 3
            error('Not sure this works for cycles of length <= 3')
        end
        
        % Cut common edge of first and last triangle
        f1 = F(cy(1), :);
        fn = F(cy(end), :);

        common_edge = intersect(f1, fn);
        
        % plot the cycle before flattening
        if debug
            figure
            mesh_cy = Mesh(V, F(cy, :));
            mesh_cy.draw('FaceAlpha', 0.8);
            p1 = V(common_edge(1), :);
            p2 = V(common_edge(2), :);
            hold on
            arrow(p1, p2, 'color', 'r')
            hold off
        end
        
        if length(common_edge) == 2
            vi = common_edge(1);
            vj = common_edge(2);

            V = [V; V(vi, :); V(vj, :)];
            vi_new = nv + 1;
            vj_new = nv + 2;
            nv = nv + 2;

            fn(fn == vi) = vi_new;
            fn(fn == vj) = vj_new;
            F(cy(end), :) = fn;
            
            % This is not very robust
            k = 1;
            fn_prev = F(cy(end-k), :);
            while ~isempty(intersect(fn_prev, [vi, vj]))
                fn_prev(fn_prev == vi) = vi_new;
                fn_prev(fn_prev == vj) = vj_new;
                F(cy(end-k), :) = fn_prev;
                k = k + 1;
                fn_prev = F(cy(end-k), :);
            end

            v1 = vi;
            v2 = vj;
            vn = vi_new;
            vn_prev = vj_new;
        elseif length(common_edge) > 2
            error('?')
        end

        % precalcuate stuff

        Vi = V(F(cy, 1), :);
        Vj = V(F(cy, 2), :);
        Vk = V(F(cy, 3), :);

        Vij = Vj - Vi;
        Vik = Vk - Vi;
        Vjk = Vk - Vj;

        Vij_norm = row_norm(Vij);
        Vik_norm = row_norm(Vik);
        Vjk_norm = row_norm(Vjk);
        norm_ratios = Vik_norm ./ Vij_norm;

        angles1 = vec_vec_angle(Vij, Vik);
        angles2 = vec_vec_angle(-Vij, Vjk);
        angles3 = vec_vec_angle(-Vik, -Vjk);

        % Build and solve system

        A = sparse(2*cy_len+4, 2*nv);
        for idx = 1:cy_len
            fid = cy(idx);
            xi = F(fid, 1);
            xj = F(fid, 2);
            xk = F(fid, 3);
            yi = xi + nv;
            yj = xj + nv;
            yk = xk + nv;

            theta = angles1(idx);
            n = norm_ratios(idx);

            A(idx, xi) = -1 + n*cos(theta);
            A(idx, xj) = -n*cos(theta);
            A(idx, xk) = 1;
            A(idx, yi) = -n*sin(theta);
            A(idx, yj) = n*sin(theta);
            A(idx, yk) = 0;

            idx2 = idx + cy_len;
            A(idx2, xi) = n*sin(theta);
            A(idx2, xj) = -n*sin(theta);
            A(idx2, xk) = 0;
            A(idx2, yi) = -1+n*cos(theta);
            A(idx2, yj) = -n*cos(theta);
            A(idx2, yk) = 1;
        end

        b = zeros(size(A, 1), 1);

        f = F(cy(1), :);
        xi = f(1);
        xj = f(2);
        yi = f(1) + nv;
        yj = f(2) + nv;

        rows = 2*cy_len;
        A(rows+1, xi) = 1;
        b(rows+1) = 0;
        A(rows+2, yi) = 1;
        b(rows+2) = 0;
        A(rows+3, xj) = 1;
        b(rows+3) = Vij_norm(1);
        A(rows+4, yj) = 1;
        b(rows+4) = 0;

        u = A \ b;
        check_norm('A*u', 'b');

        % Extract flat mesh

        Ff = F(cy, :);
        Vf = zeros(size(V));
        Vf(:, 1) = u(1:nv);
        Vf(:, 2) = u(nv+1:end);
        mesh_flat = Mesh(Vf, Ff);
        flat_meshes{end+1} = mesh_flat;

        % Plot

        if debug
            figure()
            subplot(121)
            hold on
            mesh.draw('FaceAlpha', 0.8)
            mesh.plotH(1);
            mesh.drawLabels('FontSize', 6);
            hold off
            subplot(122)
            mesh_flat.draw()
            hold on
            mesh_flat.drawLabels('FontSize', 6);
            hold off
        end

        % Check that the angles and edge lengths are the same
        
        if debug
            Vfi = Vf(F(cy, 1), :);
            Vfj = Vf(F(cy, 2), :);
            Vfk = Vf(F(cy, 3), :);

            Vfij = Vfj - Vfi;
            Vfik = Vfk - Vfi; 
            Vfjk = Vfk - Vfj;

            Vfij_norm = row_norm(Vfij);
            Vfik_norm = row_norm(Vfik);
            Vfjk_norm = row_norm(Vfjk);
            normsf = Vfik_norm ./ Vfij_norm;

            anglesf = vec_vec_angle(Vfij, Vfik);
            anglesf2 = vec_vec_angle(-Vfij, Vfjk);
            anglesf3 = vec_vec_angle(-Vfik, -Vfjk);

            assert(abs(Vfij_norm - Vij_norm) < 1e-10)
            assert(abs(Vfik_norm - Vik_norm) < 1e-10)
            assert(abs(angles1 - anglesf) < 1e-10)
            assert(abs(angles2 - anglesf2) < 1e-10)
            assert(abs(angles3 - anglesf3) < 1e-10)
        end

        % Calculate angle defect
        e1 = Vf(v2, :) - Vf(v1, :);
        en = Vf(vn, :) - Vf(vn_prev, :);
        z = cross(en, e1);
        z = sign(z(3));
        defects(i) = z*vec_vec_angle(e1, en);
        %Ad = generator_angle_defects(mesh);
        %defects(i) = vec_vec_angle(e1, en);
        
    end


end

