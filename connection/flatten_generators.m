function [defects, flat_meshes] = flatten_generators(mesh, debug)
    %function [defects, flat_meshes] = flatten_generators(mesh, debug)
    
    if nargin < 2
        debug = false;
    end

    flat_meshes = {};
    cycles = mesh.generator_cycles;
    ng2 = numel(cycles);
    defects = zeros(ng2, 1);
    
    for i = 1:ng2
        V = mesh.V; F = mesh.F;
        nv = mesh.nV; nf = mesh.nF;
        cy = cycles{i};
        cy = cy(1:end-1);
        cy_len = length(cy);
        if cy_len <= 3
            error('Not sure this works for cycles of length <= 3')
        end
        
        % Find common edge between first and last face
        f1 = F(cy(1), :);
        fn = F(cy(end), :);
        common_edge = intersect(f1, fn);
%         if length(common_edge) == 2
%             v1 = common_edge(1);
%             v2 = common_edge(2);
%             if ~mesh.directed_edge_is_in_face(cy(1), v1, v2)
%                 tmp = v1;
%                 v1 = v2;
%                 v2 = tmp;
%             end
%             
%             % Separate common edge
%             V = [V; V(v1, :); V(v2, :)];
%             v1_tilde = nv + 1;
%             v2_tilde = nv + 2;
%             nv = nv + 2;
% 
%             fn(fn == v1) = v1_tilde;
%             fn(fn == v2) = v2_tilde;
%             F(cy(end), :) = fn;
%         else
%             error('Bad cycle was given')
%         end
%         
%         % Cut the cycle
%         visited = zeros(nv, 1);
%         for j = 1:cy_len
%             f = cy(j);
%             face = F(f, :);
%             for k = 1:3
%                 v = face(k);
%                 if visited(v)
%                     nv = nv + 1;
%                     face(face == v) = nv;
%                     V = [V; V(v, :)];
%                 else
%                     visited(v) = true;
%                 end
%             end
%             F(f, :) = face;
%         end


        
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
            v1 = common_edge(1);
            v2 = common_edge(2);
            if ~mesh.directed_edge_is_in_face(cy(1), v1, v2)
                tmp = v1;
                v1 = v2;
                v2 = tmp;
            end
                
            V = [V; V(v1, :); V(v2, :)];
            v1_tilde = nv + 1;
            v2_tilde = nv + 2;
            nv = nv + 2;

            fn(fn == v1) = v1_tilde;
            fn(fn == v2) = v2_tilde;
            F(cy(end), :) = fn;
            
            % This is probably not very robust
            k = 1;
            fn_prev = F(cy(end-k), :);
            while ~isempty(intersect(fn_prev, [v1, v2]))
                fn_prev(fn_prev == v1) = v1_tilde;
                fn_prev(fn_prev == v2) = v2_tilde;
                F(cy(end-k), :) = fn_prev;
                k = k + 1;
                fn_prev = F(cy(end-k), :);
            end

%             % Cut the cycle
%             visited = zeros(nv, 1);
%             for j = 1:cy_len
%                 f = cy(j);
%                 face = F(f, :);
%                 for k = 1:3
%                     v = face(k);
%                     if visited(v)
%                         nv = nv + 1;
%                         face(face == v) = nv;
%                         V = [V; V(v, :)];
%                     else
%                         visited(v) = true;
%                     end
%                 end
%                 F(f, :) = face;
%             end

            %v1 = v1;
            %v2 = v2;
            %%vn = vi_new;
            %%vn_prev = vj_new;
            %v1_tilde = v1_tilde;
            %v2_tilde = v2_tilde;
        else
            error('Bad cycle was given')
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

%         lij = row_norm(Vij);
%         lik = row_norm(Vik);
%         ljk = row_norm(Vjk);
%         as = row_norm_square(Vij);
%         bs = row_norm_square(Vik);
%         cs = row_norm_square(Vjk);
%         
%         cos1 = 0.5 * (as + bs - cs) ./ (lij.*lik);
%         s = 0.5*(lij + lik + ljk);
%         area = s .* (s - lij) .* (s - lik) .* (s - ljk);
%         area(area<0) = 0;
%         area = sqrt(area);
%         sin1 = 2*area./(lij.*lik);
%         norm_ratios = lik ./ lij;
        
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

%             %theta = angles1(idx);
%             n = norm_ratios(idx);
% 
%             A(idx, xi) = -1 + n*cos1(idx);
%             A(idx, xj) = -n*cos1(idx);
%             A(idx, xk) = 1;
%             A(idx, yi) = -n*sin1(idx);
%             A(idx, yj) = n*sin1(idx);
%             A(idx, yk) = 0;
% 
%             idx2 = idx + cy_len;
%             A(idx2, xi) = n*sin1(idx);
%             A(idx2, xj) = -n*sin1(idx);
%             A(idx2, xk) = 0;
%             A(idx2, yi) = -1+n*cos1(idx);
%             A(idx2, yj) = -n*cos1(idx);
%             A(idx2, yk) = 1;
        end

        rhs = zeros(size(A, 1), 1);

        %f = F(cy(1), :);
        %xi = f(1);
        %xj = f(2);
        %yi = f(1) + nv;
        %yj = f(2) + nv;
        xi = v1;
        xj = v2;
        yi = xi + nv;
        yj = xj + nv;

        rows = 2*cy_len;
        A(rows+1, xi) = 1;
        rhs(rows+1) = 0;
        A(rows+2, yi) = 1;
        rhs(rows+2) = 0;
        A(rows+3, xj) = 1;
        rhs(rows+3) = Vij_norm(1);
        A(rows+4, yj) = 1;
        rhs(rows+4) = 0;

        u = A \ rhs;
        check_norm('A*u', 'rhs');

        % Extract flat mesh

        Ff = F(cy, :);
        Vf = zeros(size(V));
        Vf(:, 1) = u(1:nv);
        Vf(:, 2) = u(nv+1:end);
        mesh_flat = Mesh(Vf, Ff);
        flat_meshes{end+1} = mesh_flat;

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

            %assert(abs(Vfij_norm - Vij_norm) < 1e-10)
            %assert(abs(Vfik_norm - Vik_norm) < 1e-10)
            assert(abs(angles1 - anglesf) < 1e-10)
            assert(abs(angles2 - anglesf2) < 1e-10)
            assert(abs(angles3 - anglesf3) < 1e-10)
        end

        % Calculate angle defect
        e1 = Vf(v2, :) - Vf(v1, :);
        %en = Vf(vn_prev, :) - Vf(vn, :);
        %en = Vf(v2_tilde, :) - Vf(v1_tilde, :);
        en = Vf(v2_tilde, :) - Vf(v1_tilde, :);
        z = cross(en, e1);
        z = sign(z(3));
        defects(i) = z*vec_vec_angle(e1, en);
        %Ad = generator_angle_defects(mesh);
        %defects(i) = vec_vec_angle(e1, en);
        
        % Plot

        if debug
            figure()
            %subplot(121)
            %hold on
            %mesh.draw('FaceAlpha', 0.8)
            %mesh.plotH(1);
            %mesh.drawLabels('FontSize', 6);
            %hold off
            %subplot(122)
            p1 = Vf(v1, :)+[0,0.01,0];
            p2 = Vf(v2, :)+[0,0.01,0];
            p1_tilde = Vf(v1_tilde, :)+[0,0.01,0];
            p2_tilde = Vf(v2_tilde, :)+[0,0.01,0];
            mesh_flat.draw()
            title([num2str(defects(i)*180/pi), ' degrees'])
            hold on
            %mesh_flat.drawLabels('FontSize', 8);
            text(p1(1), p1(2), 'v1', 'color', 'r')
            text(p2(1), p2(2), 'v2', 'color', 'r')
            text(p1_tilde(1), p1_tilde(2), '\tilde v1', 'color', 'r')
            text(p2_tilde(1), p2_tilde(2), '\tilde v2', 'color', 'r')
            hold off
        end
        
    end


end

