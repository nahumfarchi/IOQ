classdef NRosy < handle
    %NROSY Create a N-RoSy directional field on the given mesh with the
    %given constraints.
    %
    % Some of the code here is based on libigl: 
    %    http://libigl.github.io/libigl/
    %
    % Example:
    %   m = Mesh();
    %   m.loadTM(path_to_off_file);
    %   N = 4;
    %   constrained_faces = [1];
    %   constraint_vectors = [1, 1, 1];
    %   solver = NRows(m, N, constrained_faces, constraint_vectors);
    %   solver.solve_with_direct_rounding();
    
    properties
        N                   % The number of vectors in the directional field.
        mesh
        local_frames        % 2*nF x 3 matrix of local frames. First nF 
                            % rows are for e1, and the last nF rows are for e2.
                            
                            % (B is the number of constraints)
        %constrained_faces   % B x 1 matrix with the indices of the constrained faces.
        %constraint_vectors  % B x 3 matrix with the constraint vectors.
        face_is_constrained % B x 1 boolean vector.
        period_is_fixed     % nE x 1 boolean vector.
        
        k                   % nE x 1 Angle difference of any two local frames.
        A                   % MIQ system of equalities Ax=b.
        b
        thetas              % nF x 1 Representative angle of the directional field
                            % each face (with respect to the local frame).
        periods             % nE x 1 
        
        n_theta             % Number of free angle variables.
        n_period            % Number of free face variables.
        n_variables         % Total number of free variables.
        
        theta_inds          % Maps the variables back into the right place in
                            % theta. I.e,
                            %   theta(theta_inds(i)) = x(i);
        period_inds 
        
        field               % N*nF x 3 Directional field. First N rows are for
                            % the first vector in each face, second N rows
                            % are for the second vector in each face, etc.
    end
    
    methods
        function obj = NRosy(mesh, N, constrained_faces, constraint_vectors)
            if N ~= 4
                warning('N != 4 is not well tested!');
            end
            obj.N = N;
            obj.mesh = mesh;
            %obj.constrained_faces = constrained_faces;
            %obj.constraint_vectors = constraint_vectors;
            obj.local_frames = create_local_frames(mesh);
            obj.compute_frame_diffs();
            
            % constant frames
            %obj.local_frames = [repmat([1,0,0], mesh.nF, 1); repmat([0,1,0], mesh.nF, 1)];
            %obj.k = zeros(mesh.nE, 1);
            
            obj.period_is_fixed = zeros(obj.mesh.nE, 1);
            obj.periods = zeros(obj.mesh.nE, 1);
            
            if nargin > 3
                obj.create_NRosy_system(constrained_faces, constraint_vectors, true);
            end
        end
        
        function compute_frame_diffs(self, debug)
            %COMPUTE_FRAME_DIFFS Compute the angles between local frames of adjacent
            %faces in the given mesh.

            if nargin < 2
                debug = false;
            end
            
            m = self.mesh;

            EV = m.EVAdj;
            EF = m.EFAdj;
            face_normals = m.FNormals;
            V = m.V;
            F = m.F;
            self.k = zeros(m.nE, 1);
            assert(m.nE == size(EF, 1))
            assert(m.nE == size(EV, 1))

            for eid = 1:m.nE
                if m.isBoundaryEdge(eid)
                    continue;
                end

                fid0 = EF(eid, 1);
                fid1 = EF(eid, 2);

                N0 = face_normals(fid0, :);
                %N1 = face_normals(fid1, :);

                v1 = EV(eid, 1);
                v2 = EV(eid, 2);

                %common_edge = normalize_rows(V(EV(eid, 2),:) - V(EV(eid, 1),:));
                common_edge = normalize_rows(V(v2, :) - ...
                                             V(v1, :));

                % Transform the first triangle so that:
                %   x axis is the common edge
                %   z axis is N0
                % the triangle now sits in the XY plane.
                o = V(EV(eid, 1),:);
                tmp = -cross(N0, common_edge, 2);
                P = [common_edge; tmp; N0];

                V0 = [V(F(fid0,1),:) - o; ...
                      V(F(fid0,2),:) - o; ...
                      V(F(fid0,3),:) - o];

                % Apply the same transformation to the second triangle.
                V1 = [V(F(fid1,1),:) - o; ...
                      V(F(fid1,2),:) - o; ...
                      V(F(fid1,3),:) - o];

                V0 = (P * V0')';
                V1 = (P * V1')';

                % Make sure all z-coordinates are zero.
                assert(V0(1,3) < 1e-10)
                assert(V0(2,3) < 1e-10)
                assert(V0(3,3) < 1e-10)

                if debug
                    figure;
                    subplot(211);
                    self.debug_plot(V0, V1, P);
                    title('Before rotation');
                end

                % Rotate the second triangle so that it sits in XY plane as well.
                % We have to find the vertex that is not on the common edge.
                other_vert_row = -1;
                for i = 1:3
                    if F(fid1, i) ~= EV(eid, 1) && F(fid1, i) ~= EV(eid, 2)
                        other_vert_row = i;
                    end
                end
                alpha = -atan2(V1(other_vert_row, 3), V1(other_vert_row, 2));
                R = [1, 0, 0; ...
                     0, cos(alpha), -sin(alpha); ...
                     0, sin(alpha), cos(alpha)];
                V1 = (R*V1')';

                if debug
                    subplot(212);
                    self.debug_plot(V0, V1, P);
                    title('Rotated');
                end

                % Make sure V1 now sits in the XY plane
                assert(V1(1,3) < 1e-10);
                assert(V1(2,3) < 1e-10);
                assert(V1(3,3) < 1e-10);

                % Calculate the angle between the two local frames.
                % Reminder: for each frame, e0 is the first edge 
                % of the triangle
                f0_e0 = normalize_rows(V0(2, :) - V0(1, :));
                f1_e0 = normalize_rows(V1(2, :) - V1(1, :));
                angle_between_frames = atan2(f1_e0(2), f1_e0(1)) - atan2(f0_e0(2), f0_e0(1));

                % Make sure that rotation is correct.
                R = [cos(angle_between_frames), -sin(angle_between_frames); ...
                     sin(angle_between_frames), cos(angle_between_frames)];
                f0_e0_rot = (R*f0_e0(1:2)')';
                assert(norm(f0_e0_rot - f1_e0(1:2)) < 1e-10)

                % It's all good.
                self.k(eid) = angle_between_frames;
            end
        end

        function debug_plot(self, V0, V1, P)
            % Draw the two triangles
            F0tmp = [1, 2, 3];
            F1tmp = [1, 2, 3];
            patch('faces', F0tmp, 'vertices', V0, 'FaceColor', 'blue')
            hold on
            patch('faces', F1tmp, 'vertices', V1, 'FaceColor', 'blue')
            for i = 1:3
                x = V0(i, 1); y = V0(i, 2); z = V0(i, 3);
                text(x, y, z, num2str(i))
                x = V1(i, 1)+0.1; y = V1(i, 2); z = V1(i, 3);
                text(x, y, z, num2str(i), 'Color', 'red')
            end
            % Draw the XYZ axis
            quiver3(0, 0, 0, P(1,1), P(1,2), P(1,3), 'Color', 'red')
            quiver3(0, 0, 0, P(2,1), P(2,2), P(2,3), 'Color', 'green')
            quiver3(0, 0, 0, P(3,1), P(3,2), P(3,3), 'Color', 'blue')
            xlabel('x'); ylabel('y'); zlabel('z');
            hold off
        end
        
        function create_NRosy_system(self, constrained_faces, constraint_vectors, reduce_space)
            %CREATE_NROSY_SYSTEM Create mixed interger system for N-Rosy field.
            % constrained_faces  - B x 1 fids of the constrained faces.
            % constraint_vectors - B x 3 vectors for the constrained faces.
            
            if nargin < 4
                reduce_space = true;
            end
            
            m = self.mesh;

            B = length(constrained_faces);

            %self.local_frames = create_local_frames(m);

            % Project each constraint vector onto its reference frame and
            % calculate the angle.
            %self.thetas = -ones(m.nF, 1);
            self.thetas = zeros(m.nF, 1);
            self.face_is_constrained = zeros(m.nF, 1);
            for i = 1:length(constrained_faces)
                fid = constrained_faces(i);
                v = constraint_vectors(i, :);
                frame = [self.local_frames(fid, :);
                         self.local_frames(fid+m.nF, :)];
                v_proj = frame * v';
                %self.thetas(fid) = atan2(v_proj(1), v_proj(2));
                self.thetas(fid) = atan2(v_proj(2), v_proj(1));
                self.face_is_constrained(fid) = true;
            end

            EV = m.EVAdj;
            EF = m.EFAdj;
            assert(size(EV,1) == m.nE);

            % Figure out how many free variables there are.
            self.n_theta = m.nF - B;
            count = 1;
            self.theta_inds = -ones(m.nF, 1);
            for fid = 1:m.nF
                if ~self.face_is_constrained(fid)
                    self.theta_inds(fid) = count;
                    count = count + 1;
                end
            end
            assert(self.n_theta == count - 1);

            disp(['n_theta : ', num2str(self.n_theta)])

            self.period_inds = -ones(m.nE, 1);
            %self.period_is_fixed = zeros(m.niE, 1);
            %self.periods = zeros(m.niE, 1);
            if reduce_space
                [self.period_is_fixed, self.periods] = reduce_search_space(...
                    m, ...
                    self.face_is_constrained, ...
                    self.thetas, ...
                    self.k);
            end
            %else
            %    self.period_is_fixed = zeros(m.nE, 1);
            %    self.periods = zeros(m.nE, 1);
            %end

            self.n_period = length(self.period_is_fixed) - sum(self.period_is_fixed);
            for eid = 1:m.nE
                if ~self.period_is_fixed(eid)
                    self.period_inds(eid) = count;
                    count = count + 1;
                end
            end

            disp(['n_period : ', num2str(self.n_period)])

            self.n_variables = self.n_theta + self.n_period;
            assert(self.n_variables == count - 1);
            self.b = zeros(self.n_variables, 1);
            I = [];
            J = [];
            V = [];

            disp(['n_variables : ', num2str(self.n_variables)])

               
            for eid = 1:m.nE
                fid1 = EF(eid, 1);
                fid2 = EF(eid, 2);
                
                if fid1 <= 0 || fid2 <= 0 || (self.face_is_constrained(fid1) && self.face_is_constrained(fid2))
                    continue
                end
                
                if fid1 >= 1
                    theta1 = self.thetas(fid1);
                    t1_index = self.theta_inds(fid1);
                end
                if fid2 >= 1
                    theta2 = self.thetas(fid2);
                    t2_index = self.theta_inds(fid2);
                end
                
                %disp(eid-1)
             
                k12 = self.k(eid);
                p12 = self.periods(eid);
                p_index = self.period_inds(eid);
                
                isFixed1 = self.face_is_constrained(fid1);
                isFixed2 = self.face_is_constrained(fid2);
                isFixedp = self.period_is_fixed(eid);

                if ~isFixed1
                    row = t1_index;
                    if isFixed1 
                        self.b(row) = self.b(row) - 2 * theta1;
                    else
                        I(end+1) = row;
                        J(end+1) = t1_index;
                        V(end+1) = 2;
                    end
                    if isFixed2
                        self.b(row) = self.b(row) + 2 * theta2;
                    else
                        I(end+1) = row;
                        J(end+1) = t2_index;
                        V(end+1) = -2;
                    end
                    if isFixedp
                        self.b(row) = self.b(row) - (4*pi / self.N) * p12;
                    else
                        I(end+1) = row;
                        J(end+1) = p_index;
                        V(end+1) = 4*pi / self.N;
                    end
                    self.b(row) = self.b(row) - 2*k12;
                end
                
                if ~isFixed2
                    row = t2_index;
                    if isFixed1 
                        self.b(row) = self.b(row) + 2 * theta1;
                    else
                        I(end+1) = row;
                        J(end+1) = t1_index;
                        V(end+1) = -2;
                    end
                    if isFixed2
                        self.b(row) = self.b(row) - 2 * theta2;
                    else
                        I(end+1) = row;
                        J(end+1) = t2_index;
                        V(end+1) = 2;
                    end
                    if isFixedp
                        self.b(row) = self.b(row) + (4*pi / self.N) * p12;
                    else
                        I(end+1) = row;
                        J(end+1) = p_index;
                        V(end+1) = -4*pi / self.N;
                    end
                    self.b(row) = self.b(row) + 2*k12;
                end
                
                if ~isFixedp
                    row = p_index;
                    if isFixed1 
                        self.b(row) = self.b(row) - (4*pi/self.N) * theta1;
                    else
                        I(end+1) = row;
                        J(end+1) = t1_index;
                        V(end+1) = (4*pi)/self.N;
                    end
                    if isFixed2
                        self.b(row) = self.b(row) + (4*pi/self.N) * theta2;
                    else
                        I(end+1) = row;
                        J(end+1) = t2_index;
                        V(end+1) = -(4*pi)/self.N;
                    end
                    if isFixedp
                        self.b(row) = self.b(row)  - 2*((2*pi/self.N)^2) * p12;
                    else
                        I(end+1) = row;
                        J(end+1) = p_index;
                        V(end+1) = 2*((2*pi/self.N)^2);
                    end
                    self.b(row) = self.b(row) - (4 * pi)/self.N * k12;
                end
                
%                 if fid1 >= 1 && ~self.face_is_constrained(fid1)
%                     I(end+1) = t1_index;
%                     J(end+1) = t1_index;
%                     V(end+1) = 2;
% 
%                     self.b(t1_index) = self.b(t1_index) - 2*k12;
% 
%                     if self.period_is_fixed(eid)
%                         self.b(t1_index) = self.b(t1_index) - 4*pi*p12/self.N;
%                     else
%                         I(end+1) = t1_index;
%                         J(end+1) = p_index;
%                         V(end+1) = 4*pi/self.N;
%                         assert(p_index > 0);
%                     end
%                     
%                     if fid2 >= 1
%                         if self.face_is_constrained(fid2)
%                             self.b(t1_index) = self.b(t1_index) + 2*theta2;
%                         else
%                             I(end+1) = t1_index;
%                             J(end+1) = t2_index;
%                             V(end+1) = -2;
%                             assert(t2_index > 0);
%                         end
%                     end
% 
%                     assert(t1_index > 0);
%                 end
% 
%                 if fid2 >= 1 && ~self.face_is_constrained(fid2)
%                     if fid1 >= 1
%                         if self.face_is_constrained(fid1)
%                             self.b(t2_index) = self.b(t2_index) + 2*theta2;
%                         else
%                             I(end+1) = t2_index;
%                             J(end+1) = t1_index;
%                             V(end+1) = -2;
%                             assert(t1_index > 0);
%                         end
%                     end
% 
%                     self.b(t2_index) = self.b(t2_index) + 2*k12;
% 
%                     if self.period_is_fixed(eid)
%                         % TODO
%                         self.b(t2_index) = self.b(t2_index) + 4*pi*p12/self.N;
%                     else
%                         I(end+1) = t2_index;
%                         J(end+1) = p_index;
%                         % TODO Why did I put a minus here?
%                         V(end+1) = -4*pi/self.N;
%                         assert(p_index > 0);
%                     end
% 
%                     I(end+1) = t2_index;
%                     J(end+1) = t2_index;
%                     V(end+1) = 2;
%                     assert(t2_index > 0);
%                 end
% 
%                 if ~self.period_is_fixed(eid)
%                     if fid1 >= 1
%                         if self.face_is_constrained(fid1)
%                             self.b(p_index) = self.b(p_index) - 4*pi*theta1/self.N;
%                         else
%                             I(end+1) = p_index;
%                             J(end+1) = t1_index;
%                             V(end+1) = 4 * pi / self.N;
%                             assert(t1_index > 0);
%                         end
%                     end
% 
%                     self.b(p_index) = self.b(p_index) - 4*pi*k12/self.N;
% 
%                     I(end+1) = p_index;
%                     J(end+1) = p_index;
%                     V(end+1) = 8*pi^2/self.N^2;
%                     %V(end+1) = 8 * pi^2 / self.N;
%                     %V(end+1) = pi^2 / 2;
%                     assert(p_index > 0);
%                     
%                     if fid2 >= 1
%                         if self.face_is_constrained(fid2)
%                             self.b(p_index) = self.b(p_index) + 4*pi*theta2/self.N;
%                         else
%                             I(end+1) = p_index;
%                             J(end+1) = t2_index;
%                             V(end+1) = -4*pi/self.N;
%                             assert(t2_index > 0);
%                         end
%                     end
%                 end
            end

            self.A = sparse(I, J, V, self.n_variables, self.n_variables);
        end
        
        function solve_with_direct_rounding(self)
            m = self.mesh;
            x = self.A \ self.b;
            % Round all integer variables
            x(self.n_theta+1:end) = round(x(self.n_theta+1:end), 0);
            
            for fid = 1:m.nF
                if ~self.face_is_constrained(fid)
                    i = self.theta_inds(fid);
                    self.thetas(fid) = x(i);
                end
            end
            
            for eid = 1:m.nE
                if ~self.period_is_fixed(eid)
                    i = self.period_inds(eid);
                    self.periods(eid) = x(i);
                end
            end
        end
        
        function edgelist_phase_unwrapping(self, ub, lb)
            if nargin < 2, ub = 50; end
            if nargin < 3, lb = -50; end
            
            m = self.mesh;
            dual_n_nodes = m.nF; % Dual graph
            dual_nE = m.niE;
            wrapped = self.thetas;
              
            node_index = (1:dual_n_nodes)';
            P_index = (dual_n_nodes+1:dual_n_nodes+dual_nE)';
            Q_index = (dual_n_nodes+dual_nE+1:dual_n_nodes+2*dual_nE)';
            n_variables = dual_n_nodes + 2*dual_nE;

            beq = zeros(dual_nE, 1);
            I = [];
            J = [];
            V = [];
            EF = m.EFAdj;
            
            % Tag the edges of the dual graph.
            deids = sparse(dual_n_nodes, dual_n_nodes);
            count = 1;
            for eid = 1:m.nE
                di = EF(eid, 1);
                dj = EF(eid, 2);
                if di < 1 || dj < 1
                    warning('Boundary edge are not well tested yet.')
                    continue;
                end
                deids(di, dj) = count;
                count = count + 1;
            end
            
            % [n_j   n_i   P_ij   Q_ij | b                        ]
            % [-1    1     1      -1   | TODO ]
            for eid = 1:m.nE
                di = EF(eid, 1);
                dj = EF(eid, 2);
                if di < 1 || dj < 1
                    warning('Boundary edges are not well tested yet.')
                    continue;
                end
                % Dual graph
                deid = deids(di, dj);
                ni = node_index(di);
                nj = node_index(dj);
                Pij = P_index(deid);
                Qij = Q_index(deid);

                I(end+1) = deid;
                J(end+1) = ni;
                V(end+1) = 1;

                I(end+1) = deid;
                J(end+1) = nj;
                V(end+1) = -1;

                I(end+1) = deid;
                J(end+1) = Pij;
                V(end+1) = 1;

                I(end+1) = deid;
                J(end+1) = Qij;
                V(end+1) = -1;

                ti = wrapped(di);
                tj = wrapped(dj);
                kij = self.k(eid);
                period_ij = self.periods(eid);
                
                %beq(deid) = round((ti-tj+kij+period_ij)/(2*pi), 0);
                
                beq(deid) = round((tj-ti)/(2*pi), 0);
                
                %beq(deid) = round((ti-tj+kij+period_ij)/(2*pi/self.N), 0);
                
                %beq(deid) = round((ti + self.k(eid) + (2*pi/self.N)*self.periods(eid) - tj) / (2*pi/self.N), 0);
                %beq(deid) = round((ti + kij + (2*pi/self.N)*period_ij - tj) / (2*pi), 0);
                %beq(deid) = round((ti-tj+kij+(2*pi/self.N)*period_ij)/(2*pi/self.N), 0);
                %beq(deid) = round(self.N * (ti + kij + (2*pi/self.N)*period_ij - tj) / (2*pi), 0);
                
                %beq(deid) = round((ti-tj+kij) / (2*pi/self.N), 0); % (!)
            end

            assert(length(I) == 4*dual_nE)
            assert(length(J) == 4*dual_nE)
            assert(length(V) == 4*dual_nE)

            Aeq = full(sparse(I, J, V, dual_nE, n_variables));

            LB = lb*ones(n_variables, 1);
            LB(dual_n_nodes+1:end) = 0;
            UB = ub*ones(n_variables, 1);
            UB(dual_n_nodes+1:end) = 10;

            cost = ones(n_variables, 1);
            %cost(Q_index) = -1;
            %cost(1:dual_n_nodes) = 0;

            [x,fval,exitflag] = linprog(cost, [], [], Aeq, beq, LB, UB);
            %[x,fval,exitflag] = intlinprog(cost, 1:dual_n_nodes, [], [], Aeq, beq, LB, UB);
            %x1 = round(x1, 0);
            % x(1:dual_n_nodes) - n_i
            % x(dual_n_nodes+1:end) - P_ij, Q_ij
            x = round(x, 0);

            unwrapped = zeros(dual_n_nodes, 1);
            for di = 1:dual_n_nodes
               %unwrapped(di) = wrapped(di) + 2*pi*x(di)/self.N;
               unwrapped(di) = wrapped(di) + 2*pi*x(di)/self.N;
            end
            
            disp(['theta diff: ', num2str(norm(mod(self.thetas,2*pi) - mod(unwrapped,2*pi)))]);
            
            self.thetas = mod(unwrapped, 2*pi);
            
            constrained_faces = 1:m.nF;
            constraint_vectors = zeros(m.nF, 3);
            for fid = 1:m.nF
                frame = [self.local_frames(fid, :); ...
                         self.local_frames(fid+m.nF, :)];
                t = self.thetas(fid);
                vec = [cos(t), sin(t)] * frame;
                constraint_vectors(fid, :) = vec;
            end

%             for eid = 1:m.nE
%                 if ~m.isBoundaryEdge(eid)
%                     continue;
%                 end
%                 
%                 fid1 = EF(eid, 1);
%                 fid2 = EF(eid, 2);
%                 ni = x(fid1);
%                 nj = x(fid2);
%                 self.periods(eid) = self.periods(eid) + ni - nj;
%             end
%             
%             constrained_faces = [];
%             constraint_vectors = [];
            
            self.create_NRosy_system(constrained_faces, constraint_vectors, false);
            self.solve_with_direct_rounding();
        end
        
        function edgelist_phase_unwrapping2(self)
            m = self.mesh;
            nE_dual = m.niE;
            nV_dual = m.nF;
            n_rows = nE_dual;
            n_cols = 2*nE_dual + nV_dual;
            offset = 2*nE_dual;
            wrapped = self.thetas;
            
            EF = m.EFAdj;
            
            I = [];
            J = [];
            V = [];
            cost = zeros(n_cols, 1);
            lb = zeros(n_cols, 1);
            ub = zeros(n_cols, 1);
            beq = zeros(n_rows, 1);
            row = 1;
            for eid = 1:m.nE
                if m.isBoundaryEdge(eid)
                    continue;
                end
                
                fid1 = min(EF(eid, 1), EF(eid, 2));
                fid2 = max(EF(eid, 1), EF(eid, 2));
                
                % n_i
                I(end+1) = row;
                J(end+1) = offset + fid1;
                V(end+1) = 1;
                % n_j
                I(end+1) = row;
                J(end+1) = offset + fid2;
                V(end+1) = -1;
                % P_ij
                col = 2*(row-1) + 1;
                I(end+1) = row;
                J(end+1) = col;
                V(end+1) = 1;
                cost(col) = 1;
                lb(col) = 0;
                ub(col) = 10;
                % Q_ij
                col = 2*(row-1) + 2;
                I(end+1) = row;
                J(end+1) = col;
                V(end+1) = -1;
                cost(col) = 1;
                lb(col) = 0;
                ub(col) = 10;      
                % RHS
                t1 = wrapped(fid1);
                t2 = wrapped(fid2);
                kij = self.k(eid);
                period_ij = self.periods(eid);
                beq(row) = round((t1-t2+kij)/(2*pi/self.N), 0);
                
                row = row + 1;
            end
            
            cost(offset+1:end) = 0;
            lb(offset+1:end) = -50;
            ub(offset+1:end) = 50;
            
            Aeq = full(sparse(I, J, V, n_rows, n_cols));
            
            [x,fval,exitflag] = linprog(cost, [], [], Aeq, beq, lb, ub);
            x = round(x, 0);

%             unwrapped = zeros(m.nF, 1);
%             for fid = 1:m.nF
%                %unwrapped(di) = wrapped(di) + 2*pi*x(di)/self.N;
%                unwrapped(fid) = wrapped(fid) + 2*pi*x(offset+fid)/self.N;
%             end
%             
%             disp(['diff: ', num2str(norm(mod(self.thetas,2*pi) - mod(unwrapped,2*pi)))]);
%             
%             self.thetas = unwrapped;
            
            for eid = 1:m.nE
                fid1 = EF(eid, 1);
                fid2 = EF(eid, 2);
                %self.periods(eid) = self.periods(eid) + (x(offset+fid1) - x(offset+fid2));
                self.periods(eid) = -(x(offset+fid1) - x(offset+fid2));
            end          
        end
        
        function CF = create_cross_field(self)
            m = self.mesh;
            CF = zeros(m.nF, 3);
            for fid = 1:m.nF
                frame = [self.local_frames(fid, :); ...
                         self.local_frames(fid+m.nF, :)];
                for i = 0:(self.N-1)
                    t = self.thetas(fid) + i*2*pi/self.N;
                    vec = [cos(t), sin(t)] * frame;
                    CF(fid+i*m.nF, :) = vec;
                end
            end
        end
        
        function CF = get.field(self)
            CF = self.create_cross_field();
        end
        
        function draw_constraints(self, constrained_faces, constraint_vectors)
            m = self.mesh;
            P = (m.V(m.F(constrained_faces, 1), :) + ...
                 m.V(m.F(constrained_faces, 2), :) + ...
                 m.V(m.F(constrained_faces, 3), :)) / 3;
            x = P(:, 1);
            y = P(:, 2);
            z = P(:, 3);
            u = constraint_vectors(:, 1);
            v = constraint_vectors(:, 2);
            w = constraint_vectors(:, 3);
            quiver3(x, y, z, u, v, w, 'Marker', 'o', 'AutoScale', 'on', 'AutoScaleFactor', m.avg_length/3, 'color', 'm', 'LineWidth', 2);
        end
        
        function draw_field(self, scale, cf_colors)
            m = self.mesh;
            if nargin < 2
                scale = m.avg_length;
            end
            if nargin < 3
                cf_colors = {'r', 'b', 'b', 'b'};
            end
            for i = 0:self.N-1
                % First N rows are for the first vector, second N rows
                % are for the second vector, etc.
                DF = self.field(1+i*m.nF:m.nF+i*m.nF, :);
                m.drawFaceField(DF, 'color', cf_colors{i+1}, 'AutoScale', 'on', 'AutoScaleFactor', scale)
            end
        end
        
        function E = E_quadratic(self)
            m = self.mesh;
            E = 0;
            EF = m.EFAdj;
            for eid = 1:m.nE
                fidi = EF(eid, 1);
                fidj = EF(eid, 2);
                if fidi < 1 || fidj < 1
                    continue;
                end
                ti = self.thetas(fidi);
                tj = self.thetas(fidj);
                pij = self.periods(eid);
                kij = self.k(eid);
                E = E + (ti + kij + (2*pi/self.N)*pij - tj)^2;
                %fprintf('%d : %f\n', eid, (ti + kij + (2*pi/self.N)*pij - tj)^2); 
                %E = E + (-ti - kij - (2*pi/self.N)*pij + tj)^2;
            end
        end
        
        function save(self, path)
            m = self.mesh;
            DF = self.field;
            x = m.IF2V * DF(1:m.nF, :) ./ repmat(sum(m.VVAdj, 2), 1, 3);
            fid = fopen(path, 'w');
            %fprintf(fid, '%s\n', obj.header);
            %fprintf(fid, '%d %d\n', size(DF, 1), self.N);
            %fprintf(fid, '%.6g %.6g %.6g \n', DF');
            fprintf(fid, '%d %d\n', m.nV, self.N);
            fprintf(fid, '%.6g %.6g %.6g \n', x');
            fclose(fid);
        end
    end
    
end

