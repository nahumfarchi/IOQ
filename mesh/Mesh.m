classdef Mesh < matlab.mixin.Copyable
    %MESH 3D triangle mesh represented as (V,F), where V is a 
    % nv x 3 matrix of vertex positions, and F is a nf x 3 matrix of 
    % indices into V. Each row (i,j,k) of F is a triangle with vertex 
    % positions V(i,:), V(j,:), V(k,:).
    %
    % Methods:
    %         function obj = Mesh(varargin)      
    %         function emptyCache(self)      
    %         function set.V(self, V)
    %         function set.F(self, F)
    %   Saving / Loading
    %         function loadTM(self, path)
    %         function loadOFF(self, path)
    %         function loadFField(self, path)
    %         function saveFField(self, path)
    %         function saveNRosy(self, path, face_based)
    %         function saveSingInteger(self, path)
    %         function loadOBJ(self, path)
    %         function saveTM(self, path)   
    %         function saveOFF(self, path)         
    %   Geometry / Topology
    %         function AF = get.AF(self)
    %         function AV = get.AV(self)
    %         function Gv = get.Gv(self)
    %         function Gv_inv = get.Gv_inv(self)
    %         function Gf = get.Gf(self)
    %         function Gf_inv = get.Gf_inv(self)
    %
    %         function IF2V = get.IF2V(self)
    %         function IV2F = get.IV2F(self)
    %
    %         function VVAdj = get.VVAdj(self)
    %         function compute_face_topology(self)
    %         function FE = get.FEAdj(self)       
    %         function FF = get.FFAdj(self)
    %         function FFo = get.FFo(self)       
    %         function compute_edge_topology(self)
    %         function EFAdj = get.EFAdj(self)
    %         function edge_exists = edgeExists(self, eid)
    %         function EVAdj = get.EVAdj(self)
    %         function VFAdj = get.VFAdj(self)
    %         function vf_1ring = get.vf_1ring(self)
    %         function cycles = vert_cycles(self)
    %         function FNormals = get.FNormals(self)
    %         function VNormals = get.VNormals(self)
    %         function VF_angles = get.VF_angles(self)
    %         function avg = get.avg_length(self)
    %         function vert_vects = face_vect_to_vertex_vect(self, face_vects)
    %
    %         function be = isBoundaryEdge(self, eid)
    %         function be = get.boundaryEdge(self)
    %         function bv = isBoundaryVertex(self, vid)
    %         function res = oriented_edge_is_in_face(self, fid, eid)
    %         function res = directed_edge_is_in_face(self, fid, vid1, vid2)
    %
    %         function genus = get.genus(self)   
    %         function nbE = get.nbE(self)
    %         function nE = get.nE(self)
    %         function niE = get.niE(self)
    %
    %         function valence = getValence(obj)
    %         function Kg = get.gaussian_cur(self)
    %         function gen_def = get.generator_defects(self)
    %         function Mk = get.mean_cur(self)
    %
    %         function T = get.primalMST(self)
    %         function T = get.dualMST(self)
    %         function Ad = cycle_defects(self, cycles, verbose)
    %         function cycles = get.generator_cycles(self)
    %         function alpha = cycle_to_1form(self, cycle)
    %         function C = dual_cycles_to_1forms(self, cycles)
    %         function H = get.H(self)
    %
    %         function set_ffield(self, ...
    %                 degree, ...
    %                 local_frames, ...
    %                 frame_diffs, ...
    %                 ffield_angles, ...
    %                 ffield_vectors, ...
    %                 vert_sing, ...
    %                 gen_sing, ...
    %                 miq_energy, ...
    %                 connection)
    %   Plotting
    %         function colorFaces(self, faces, color)
    %         function colorEdges(self, edges, color)
    %         function A = montage(self, filename, funcs, rows, titles, varargin)
    %         function draw(obj, varargin)
    %         function drawFaceField(self, vectors, varargin)
    %         function drawField(self, vectors, varargin)
    %         function drawFNormals(self, varargin)
    %         function drawVNormals(self, varargin)
    %         function drawLabels(self, varargin)
    %         function labelFaces(self, f_labels, FORMAT, offset, varargin)
    %         function labelEdges(self, inds, labels, FORMAT, offset, varargin)
    %         function labelVertices(self, v_labels, FORMAT, offset, varargin)     
    %         function plotH(self, inds)
    %         function plotEdgePaths(self, paths)
    
    properties
        header
        file_path
        
        nV           % Number of vertices.
        nF           % Number of faces.
        nE           % Number of edges.
        nbE          % Number of boundary edges.
        niE          % Number of internal edges.
        V            % n  x 3 vertex matrix.
        F            % nF x 3 face matrix.
        
        % Mesh topology
        
        % Adjacency matrices. Note that for some matices, each entry is a 
        % boolean value signifying whether the two elements are adjacent,
        % while for others it's an id of the adjacent element. For example,
        % a row in EFAdj has two entries with the ids of the adjacent
        % faces (with a -1 if it's a boundary edge).
        VVAdj        % nv x nv boolean vertex-vertex adjacency matrix.
        EFAdj        % nE x 2  ids edge-face adjacency matrix.
        EVAdj        % nE x 2  ids edge-vertex adjacency matrix.
        VFAdj        % nV x nF boolean vertex-face adjacency matrix.
        FEAdj        % nF x 3  ids face-edge adjacency matrix.
        FFAdj        % nF x 3  ids face-face adjacency matrix.
        FFo          % nF x nF face-face adjacency matrix with the (signed)
                     %         edge id of the common edge in each entry,
                     %         where the sign signifies the orientation.
                     
        IF2V         % Interpolation matrix from F to V.
        IV2F         % Interpolation matrix from V to F.
        
        genus        % Number of "handles".
        
        vf_1ring     % Cell array of size nV with indices of adjacent faces in
                     % each cell.
        VF_angles    % Vertex-face adjacency matrix with the 1-ring angles (normalized).
        avg_length   % Average edge length.
        boundaryEdge

        primalMST
        dualMST
        generator_cycles % Cell array of size 2g with non-contractible cycles (face-id based)
        H                % E x 2g matrix with non-contractible cycles

        % Mesh geometry
        
        % Area information
        AF           % nf x 1 Face areas.
        AV           % nv x 1 Vertex areas (defined as the sum of the face 
                     % areas in the one-ring divided by 3).
        Gv           % Vertex areas.
        Gv_inv       % Vertex areas inverse.
        Gf           % Face areas replicated three times.
        Gf_inv       % Face areas inverse.

        % Curvature
        gaussian_cur % Calculated per vertex as Kg(i)=(2pi-sum(theta_j))
                     % where the theta_j's are the angles around the
                     % vertex.
        mean_cur
        generator_defects
        
        % Normals
        FNormals     % Face normals.
        VNormals     % Vertex normals.
        
        % face-based representation of a vector field
        degree          %                 The field's degree (usually 4)
        local_frames    % nf x 3          Local frame per face (i.e., two orthogonal unit vectors)
        frame_diffs     %                 (Directed) angle difference between adjacent frames (called k in MIQ)
        ffield_angles   % nf x 1          Vector field angles relative to the local frames (called theta in MIQ)
        ffield_vectors  % (nf*degree) x 3 The field's vectors in global coordinates
        n_vert_sing     %                 Number of vertex singularities (called k in TCODS)
        n_gen_sing      %                 Number of generator singularities
        n_singularities %
        vert_sing       % nv x 2          (vi, si), where vi is the vertex, and si is the singularity index.
        gen_sing        % 2*genus x 2
        miq_energy      %                 The "smoothness" of the field (see Mixed Integer Quadrangulation)
        connection      % ne x 1          See Trivial Connections on Discrete Surfaces (Crane et al)
        periods         % ne x 1          The variable p in MIQ
    end
    
    properties(Access = private)
        % Cache some of the costly (computation wise) stuff. Each getter 
        % first checks whether the field exists in the cache, only 
        % calculating it if it doesn't. The cache is emptied whenever a new 
        % mesh is loaded.
        cached_nE
        cached_nbE
        cached_niE
        cached_VVAdj
        cached_EFAdj
        cached_VFAdj
        cached_AF
        cached_AV
        cached_IF2V
        cached_IV2F
        cached_FNormals
        cached_VNormals
        cached_genus
        cached_Gv
        cached_Gv_inv
        cached_Gf
        cached_Gf_inv
        cached_VF_angles
        cached_EVAdj
        cached_FEAdj
        cached_FFAdj
        cached_boundaryEdge
        cached_H
        cached_generator_cycles
        cached_Kg
        cached_FFo
        cached_generator_defects
    end
    
    methods
        
        function obj = Mesh(varargin)
            % Initialize mesh with another mesh:
            %   Mesh(m)
            % Initialize mesh with vertices, faces, and edges:
            %   Mesh(V, F)
            obj.emptyCache();
            if nargin == 1
                if ischar(varargin{1})
                    obj.loadTM(varargin{1});
                else
                    m = varargin{1};
                    obj.file_path = m.file_path;
                    obj.V = m.V;
                    obj.F = m.F;
                    obj.nV = m.nV;
                    obj.nF = m.nF;
                    obj.nE = m.nE;
                end
            elseif nargin == 2
                V = varargin{1};
                F = varargin{2};
                obj.file_path = [];
                obj.V = V;
                obj.F = F;
                obj.nV = size(V, 1);
                obj.nF = size(F, 1);
            end
        end
        
        function emptyCache(self)
            self.cached_nE = [];
            self.cached_nbE = [];
            self.cached_niE = [];
            self.cached_VVAdj = [];
            self.cached_EFAdj = [];
            self.cached_VFAdj = [];
            self.cached_AF = [];
            self.cached_AV = [];
            self.cached_IF2V = [];
            self.cached_IV2F = [];
            self.cached_FNormals = [];
            self.cached_VNormals = [];
            self.cached_genus = [];
            self.cached_Gv = [];
            self.cached_Gv_inv = [];
            self.cached_Gf = [];
            self.cached_Gf_inv = [];
            self.cached_VF_angles = [];
            self.cached_EVAdj = [];
            self.cached_FEAdj = [];
            self.cached_FFAdj = [];
            self.cached_boundaryEdge = [];
            self.cached_H = [];
            self.cached_generator_cycles = [];
            self.cached_Kg = [];
            self.cached_FFo = [];
            self.cached_generator_defects = [];
        end
        
        function set.V(self, V)
            self.V = V;
            self.emptyCache();
        end
        
        function set.F(self, F)
            self.F = F;
            self.emptyCache();
        end
                
        %%%%%%%%%%%%%%%%%%%%
        % Loading / saving %
        %%%%%%%%%%%%%%%%%%%%
        
        function loadTM(self, path)
            [~, ~, ext] = fileparts(char(path));
            switch ext
                case '.off'
                    self.loadOFF(path);
                case '.obj'
                    self.loadOBJ(path);
                case '.ffield'
                    self.loadFField(path);
                otherwise
                    error('Unkown file extension');
            end            
            assert(size(self.V, 2) == 3, 'V does not have 3 columns!')
        end
        
        function loadOFF(self, path)
            % OFF format:
            % OFF
            % number of vertices, number of faces, number of edges
            % list of vertices X,Y,Z coordinates
            % list of faces: zero indexed
            %
            % For example:
            % OFF
            % 8 6 12
            % 1.0   0.0   1.0
            % 0.0   1.0   1.0
            % -1.0   0.0   1.0
            % 0.0  -1.0   1.0
            % 1.0   0.0  -1.0
            % 0.0   1.0  -1.0
            % -1.0   0.0  -1.0
            % 0.0  -1.0 -1.0
            % 3  0 1 2
            % 3  7 4 0 
            % 3  4 5 1 
            % 3  5 6 2 
            % 3  3 2 6  
            % 3  6 5 4
            
            self.emptyCache();
            
            self.file_path = path;
                      
            fid = fopen(path, 'r');
            self.header = fgetl(fid);
            sizes = fscanf(fid, '%d %d %d', [1, 3]); 
            self.nV = sizes(1);
            self.nF = sizes(2);
            self.cached_nE = [];
            [self.V, countV] = fscanf(fid, '%f %f %f', [3 self.nV]); 
            [self.F, countF] = fscanf(fid, '%d %f %f %f', [4 self.nF]);
            assert(countV/3 == self.nV);
            assert(countF/4 == self.nF);
            % fscanf fills obj.V and obj.F column wise so we have to transpose
            self.V = self.V';
            self.F = self.F';
            self.F = self.F(:, 2:end); % get rid of number of vertices since it's always 3
            self.F = self.F+1; % add 1 since OFF is zero-indexed
            
            fclose(fid);
        end
        
        function loadFField(self, path)
            % Load a .ffield mesh (made up format, see saveFField).
            
            self.emptyCache();
            self.file_path = path;
            fid = fopen(path, 'r');
            
            self.header = fgetl(fid);
            self.miq_energy = fscanf(fid, 'miq_energy %f\n');
            sizes = fscanf(fid, '%d %d %d\n', [1, 3]); 
            self.nV = sizes(1);
            self.nF = sizes(2);
            
            str = fgetl(fid);
            [self.V, countV] = fscanf(fid, '%f %f %f\n', [3 self.nV]); 
            assert(countV/3 == self.nV);            
            self.V = self.V';
            
            str = fgetl(fid);
            [self.F, countF] = fscanf(fid, '%d %f %f %f\n', [4 self.nF]);
            assert(countF/4 == self.nF);
            self.F = self.F';
            self.F = self.F(:, 2:end); % get rid of number of vertices since it's always 3
            self.F = self.F+1; % add 1 since OFF is zero-indexed
            
            self.degree = fscanf(fid, 'degree %d\n');
            
            str = fgetl(fid);
            [self.local_frames, n] = fscanf(fid, '%f %f %f\n', [3 2*self.nF]);
            assert(n == 3*2*self.nF);
            self.local_frames = self.local_frames';
            
            n_frame_diffs = fscanf(fid, 'frame_diffs %d\n');
            [self.frame_diffs, n] = fscanf(fid, '%f\n', [1 n_frame_diffs]);
            self.frame_diffs = self.frame_diffs(:);
            
            str = fgetl(fid);
            [self.ffield_angles, n] = fscanf(fid, '%f\n', [1 self.nF]);
            self.ffield_angles = self.ffield_angles(:);
            
            str = fgetl(fid);
            [self.ffield_vectors, n] = fscanf(fid, '%f %f %f\n', [3 self.degree*self.nF]);
            self.ffield_vectors = self.ffield_vectors';
            
            %self.n_vert_sing = fscanf(fid, 'singularities %d\n');
            %[self.vert_sing, n] = fscanf(fid, '%d %f\n', [2 self.n_vert_sing]);
            %self.vert_sing = self.vert_sing';

            %self.n_gen_sing = 0;
            %self.gen_sing = [];

            self.n_vert_sing = fscanf(fid, 'vert_sing %d\n');
            [self.vert_sing, n] = fscanf(fid, '%d %f\n', [2 self.n_vert_sing]);
            self.vert_sing = self.vert_sing';
            
            self.n_gen_sing = fscanf(fid, 'gen_sing %d\n');
            [self.gen_sing, n] = fscanf(fid, '%d %f\n', [2 self.n_gen_sing]);
            self.gen_sing = self.gen_sing';
            
            self.n_singularities = self.n_vert_sing;
            
            fclose(fid);
        end
        
        function saveFField(self, path)
            % Made up format that's similar to .off, but also saves the 
            % direction field. Note that faces are zero indexed!
            %
            % FFIELD
            % n_verts, n_faces, n_edges
            % vertices
            % vx, vy, vz
            % ...
            % faces # zero indexed
            % i1, i2, i3 # nf rows
            % ...
            % degree N
            % frames
            % k1x, k1y, k1z # nf rows
            % ...
            % k2x, k2y, k2z # nf rows
            % ...
            % frame_diffs n
            % r
            % ...
            % ffield_angles 
            % theta1
            % ...
            % ffield_vectors
            % v1x, v1y, v1z
            % ...
            % v2x, v2y, v2z
            % ... 
            % ...
            % vert_sing n_sing
            % ki, k
            % ...
            % miq_energy E
            
            fid = fopen(path, 'w+');
            %fprintf(fid, '%s\n', obj.header);
            fprintf(fid, 'FFIELD\n');
            fprintf(fid, 'miq_energy %.10g\n', self.miq_energy);
            fprintf(fid, '%d %d %d\n', self.nV, self.nF, 0);
            
            fprintf(fid, 'vertices\n');
            fprintf(fid, '%.10g %.10g %.10g \n', self.V');
            
            fprintf(fid, 'faces\n');
            fprintf(fid, '3 %d %d %d\n', (self.F-1)');
            
            fprintf(fid, 'degree %d\n', self.degree);
            
            fprintf(fid, 'frames\n');
            fprintf(fid, '%.10g %.10g %.10g\n', self.local_frames');
            
            fprintf(fid, 'frame_diffs %d\n', numel(self.frame_diffs));
            fprintf(fid, '%.10g\n', self.frame_diffs');
            
            fprintf(fid, 'ffield_angles\n');
            fprintf(fid, '%.10g\n', self.ffield_angles');
            
            fprintf(fid, 'ffield_vectors\n');
            fprintf(fid, '%.10g %.10g %.10g\n', self.ffield_vectors');
            
            fprintf(fid, 'vert_sing %d\n', self.n_vert_sing);
            fprintf(fid, '%d %.10g\n', self.vert_sing');
            
            fprintf(fid, 'gen_sing %d\n', self.n_gen_sing);
            fprintf(fid, '%d %.10g\n', self.gen_sing');
            
            fclose(fid);
        end
        
        function saveNRosy(self, path, face_based)
            % Save NRosy field in the following format:
            %
            % n_faces field_degree
            % vx, vy, vz # (n_verts rows)
            % ...
            
            fid = fopen(path, 'w+');
            
            if face_based
                n_rows = self.nF;
                X = self.ffield_vectors(1:n_rows, :);
            else
                n_rows = self.nV;
                X = self.IF2V * self.ffield_vectors(1:self.nF, :);
            end
            
            assert(size(X, 1) == n_rows, 'Field should have nV rows!');
            assert(size(X, 2) == 3, 'Filed should have 3 columns!');
            
            fprintf(fid, '%d %d\n', n_rows, self.degree);
            fprintf(fid, '%.10g %.10g %.10g \n', X');
            
            fclose(fid);
        end

        function saveSingInteger(self, path)
            fid = fopen(path, 'w+');
            
            fprintf(fid, 'degree %d\n', self.degree);

            fprintf(fid, 'vert_sing %d\n', self.n_vert_sing);
            fprintf(fid, '%d %d\n', self.vert_sing');
            
            fprintf(fid, 'gen_sing %d\n', self.n_gen_sing);
            fprintf(fid, '%d %d\n', self.gen_sing');
        end
               
        function loadOBJ(self, path)
            self.emptyCache();
            self.file_path = path;
            
            obj = read_obj(path);
            self.V = obj.v(:, 1:3);
            self.F = obj.f.v;
            self.nV = size(self.V, 1);
            self.nF = size(self.F, 1);
        end
        
        function saveTM(self, path)
            [~, ~, ext] = fileparts(char(path));
            if strcmp(ext, '.off')
                self.saveOFF(path);
            elseif strcmp(ext, '.obj')
                error('Saving OBJ format is not implemented')
            elseif strcmp(ext, '.ffield')
                self.saveFField(path);
            else
                errror('Unkown file extension');
            end
        end
        
        function saveOFF(self, path)
            fid = fopen(path, 'w+');
            fprintf(fid, 'OFF\n');
            fprintf(fid, '%d %d %d\n', self.nV, self.nF, 0);
            fprintf(fid, '%.10g %.10g %.10g \n', self.V');
            fprintf(fid, '3 %d %d %d\n', (self.F-1)');
            fclose(fid);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%
        % Geometry / topology %
        %%%%%%%%%%%%%%%%%%%%%%%
        
        % Areas
        
        function AF = get.AF(self)
            % Face area
            if isempty(self.cached_AF)
                V1 = self.V(self.F(:,2),:) - self.V(self.F(:,1),:);
                V2 = self.V(self.F(:,3),:) - self.V(self.F(:,1),:);
                self.cached_AF = 0.5*sqrt(sum(cross(V1, V2, 2).^2, 2));
            end
            AF = self.cached_AF;
        end
        
        function AV = get.AV(self)
            % Return vertex areas. Default is the sum of face areas in the
            % one-ring divided by three.
            if isempty(self.cached_AV)
                I = [self.F(:,1); self.F(:,2); self.F(:,3)];
                J = [1:self.nF, 1:self.nF, 1:self.nF];
                vals = [self.AF; self.AF; self.AF];
                self.cached_AV = sparse(I, J, vals) * ones(self.nF,1) ./ 3;
            end
            AV = self.cached_AV;
        end
        
        function Gv = get.Gv(self)
            % Diagonal matrix with vertex areas.
            if isempty(self.cached_Gv)
                self.cached_Gv = sparse(1:self.nV, 1:self.nV, self.AV, self.nV, self.nV);
            end
            Gv = self.cached_Gv;
        end
        
        function Gv_inv = get.Gv_inv(self)
            % Diagonal matrix with 1/vertex area
            if isempty(self.cached_Gv_inv)
                self.cached_Gv_inv = sparse(1:self.nV, 1:self.nV, 1./self.AV, self.nV, self.nV);
            end
            Gv_inv = self.cached_Gv_inv;
        end
        
        function Gf = get.Gf(self)
            % Diagonal matrix with face areas
            if isempty(self.cached_Gf)
                self.cached_Gf = sparse(1:3*self.nF, 1:3*self.nF, repmat(self.AF, 3, 1), 3*self.nF, 3*self.nF);
            end
            Gf = self.cached_Gf;
        end
        
        function Gf_inv = get.Gf_inv(self)
            % Diagonal matrix with 1/face areas
            if isempty(self.cached_Gf_inv)
                self.cached_Gf_inv = sparse(1:3*self.nF, 1:3*self.nF, repmat(1./self.AF, 3, 1), 3*self.nF, 3*self.nF);
            end
            Gf_inv = self.cached_Gf_inv;
        end
        
        % Interpolation
        
        function IF2V = get.IF2V(self)
            % Interpolation from face to vertex
            if isempty(self.cached_IF2V)
                I = [self.F(:,1); self.F(:,2); self.F(:,3)];
                J = [1:self.nF, 1:self.nF, 1:self.nF];
                vals = [...
                    self.AF ./ (3*self.AV(I(1:self.nF))); ...
                    self.AF ./ (3*self.AV(I(self.nF+1:2*self.nF))); ...
                    self.AF ./ (3*self.AV(I(2*self.nF+1:3*self.nF)))];
                    self.cached_IF2V = sparse(I, ...
                        J, ...
                        vals);
            end
            IF2V = self.cached_IF2V;
        end
        
        function IV2F = get.IV2F(self)
            % Interpolation from vertex to face
            if isempty(self.cached_IV2F)
                self.cached_IV2F = ...
                    spdiags(1./self.AF, 0, self.nF,self.nF) * ...
                            self.IF2V' * ...
                            spdiags(self.AV, 0, self.nV, self.nV);
            end
            IV2F = self.cached_IV2F;
        end
        
        function VVAdj = get.VVAdj(self)
            if isempty(self.cached_VVAdj)
                I1 = self.F(:,1);
                I2 = self.F(:,2);
                I3 = self.F(:,3);
                J1 = self.F(:,2);
                J2 = self.F(:,3);
                J3 = self.F(:,1);
                S = ones(3*self.nF, 1);
                self.cached_VVAdj = sparse(...
                    [I1; I2; I3], ...
                    [J1; J2; J3], ...
                    S, ...
                    self.nV, ...
                    self.nV);
                % TODO
                self.cached_VVAdj = self.cached_VVAdj | self.cached_VVAdj';
            end
            VVAdj = self.cached_VVAdj;
        end
        
        % Adjacency
        
        function compute_face_topology(self)
            FE = zeros(self.nF, 3);
            FF = zeros(self.nF, 3);
            for eid = 1:self.nE
                fid1 = self.EFAdj(eid, 1);
                fid2 = self.EFAdj(eid, 2);
                
                if fid1 >= 1
                    k = 1;
                    while FE(fid1, k) ~= 0
                        k = k + 1;
                    end
                    assert(k <= 3);
                    FE(fid1, k) = eid;
                    FF(fid1, k) = fid2;
                end
                
                if fid2 >= 1
                    k = 1;
                    while FE(fid2, k) ~= 0
                        k = k + 1;
                    end
                    assert(k <= 3);
                    FE(fid2, k) = eid;
                    FF(fid2, k) = fid1;
                end
            end
            
            self.cached_FEAdj = FE;
            self.cached_FFAdj = FF;
        end

        function compute_edge_topology(self)
            % nE x nF edge-face sparse adjacency matrix.
            % EF(i,j) == 1 iff edge i is adjacent to face j.
            I1 = (min(self.F(:,1), self.F(:,2))-1) * self.nV + max(self.F(:,1), self.F(:,2));
            I2 = (min(self.F(:,2), self.F(:,3))-1) * self.nV + max(self.F(:,2), self.F(:,3));
            I3 = (min(self.F(:,1), self.F(:,3))-1) * self.nV + max(self.F(:,1), self.F(:,3));
            J = (1:self.nF)';
            EF = sparse(...
                [I1; I2; I3], ...
                [J; J; J], ...
                ones(3*self.nF, 1), ...
                self.nV*self.nV, ...
                self.nF);

            % nE x 2 edge-vertex matrix. 
            EV = sparse(self.nV^2, 2);
            [edge_ids, ~] = find(EF > 0);
            vids_1 = floor((edge_ids-1) ./ self.nV) + 1;
            vids_2 = mod(edge_ids-1, self.nV)+1;
            EV(edge_ids, 1) = vids_1;
            EV(edge_ids, 2) = vids_2;

            % Throw away edges (i.e., rows) that don't exist in the mesh.
            [I, ~] = find(sum(EF, 2) > 0);
            EF = EF(I, :);
            EV = EV(I, :);
            
            self.cached_EVAdj = EV;
            
            % Set non-zero elements in EF to be equal to the column number 
            % that they are in.
            [I, J] = find(EF ~= 0);
            [I, ord] = sort(I);
            J = J(ord);
            ind = sub2ind(size(EF), I, J);
            EF(ind) = J;
            
            % Find first non-zero in each row of EF and place it in the
            % first column.
            [I, J] = find(EF ~= 0);
            firstIndex = accumarray(I, J, [size(EF,1), 1], @min);
            self.cached_EFAdj = zeros(size(EF, 1), 2);
            ind=sub2ind(size(EF), (1:size(EF,1))', firstIndex);
            self.cached_EFAdj(:, 1) = EF(ind);
            % Find second non-zero in each row of EF and place it in the
            % second column.
            EF(ind) = 0;
            [I, J] = find(EF ~= 0);
            firstIndex = accumarray(I, J, [size(EF,1), 1], @min, 1);
            ind = sub2ind(size(EF), (1:size(EF,1))', firstIndex);
            self.cached_EFAdj(:, 2) = EF(ind);
            % For boundary edges, zero out the second column
            [ind, ~] = find(sum(EF, 2) == 0);
            self.cached_EFAdj(ind, 2) = 0;
            
            self.cached_nbE = length(ind);
            self.cached_nE = size(EF, 1);
            self.cached_niE = self.cached_nE - self.cached_nbE;
            
            % Say that the common edge is (v1, v2) where v1 < v2.
            % Then in one face the edge will apepar as (v1, v2) and in the
            % other edge it will appear as (v2, v1). I want the face that
            % has (v1, v2) to be first.
            for eid = 1:self.nE
                fid0 = self.cached_EFAdj(eid, 1);
                v1 = EV(eid, 1);
                v2 = EV(eid, 2);
                swap = true;
                for i = 1:3
                    v1tmp = self.F(fid0, i);
                    i2 = i + 1;
                    if i2 > 3
                        i2 = 1;
                    end
                    v2tmp = self.F(fid0, i2);
                    if v1 == v1tmp && v2 == v2tmp
                        swap = false;
                        break
                    end
                end
                if swap
                    tmp = self.cached_EFAdj(eid, 1);
                    self.cached_EFAdj(eid, 1) = self.cached_EFAdj(eid, 2);
                    self.cached_EFAdj(eid, 2) = tmp;
                end
            end 
            
            %self.cached_EFAdj = [self.cached_EFAdj(:, 2), self.cached_EFAdj(:, 1)];
        end
        
        function FE = get.FEAdj(self)
            if isempty(self.cached_FEAdj)
                self.compute_face_topology();
            end
            FE = self.cached_FEAdj;
        end
        
        function FF = get.FFAdj(self)
            if isempty(self.cached_FFAdj)
                self.compute_face_topology();
            end
            FF = self.cached_FFAdj;
        end
               
        function FFo = get.FFo(self)
            if ~isempty(self.cached_FFo)
                self.FFo = self.cached_FFo;
                return;
            end
            
            % Map adjacent faces to the common edge, with the sign 
            % signifying the orientation.
            EF = self.EFAdj;
            FFo = sparse(self.nF, self.nF);
            for eid = 1:self.nE
                % EF was constructed so that the face to the left is
                % always the first face (wrt the edge orientation).
                f1 = EF(eid, 1);
                f2 = EF(eid, 2);
                FFo(f1, f2) = -eid;
                FFo(f2, f1) = eid;
            end
        end
        
        function EFAdj = get.EFAdj(self)
            if isempty(self.cached_EFAdj) || isempty(self.cached_EVAdj)
                self.compute_edge_topology();
            end
            EFAdj = self.cached_EFAdj;
        end
        
        function EVAdj = get.EVAdj(self)
            if isempty(self.cached_EFAdj) || isempty(self.cached_EVAdj)
                self.compute_edge_topology();
            end
            EVAdj = self.cached_EVAdj;
        end
        
        function VFAdj = get.VFAdj(self)
            if isempty(self.cached_VFAdj)
                I1 = self.F(:, 1);
                I2 = self.F(:, 2);
                I3 = self.F(:, 3);
                J = (1:self.nF)';
                S = ones(3*self.nF, 1);
                self.cached_VFAdj = sparse(...
                    [I1; I2; I3], ...
                    [J; J; J], ...
                    S, ...
                    self.nV, ...
                    self.nF);
            end
            VFAdj = self.cached_VFAdj;
        end

        function vf_1ring = get.vf_1ring(self)
            vf_1ring = cell(self.nV, 1);
            for i = 1:self.nF
                vf_1ring{self.F(i, 1)}(end+1) = i;
                vf_1ring{self.F(i, 2)}(end+1) = i;
                vf_1ring{self.F(i, 3)}(end+1) = i;
            end
        end
        
        function be = isBoundaryEdge(self, eid)
            if self.EFAdj(eid, 1) < 1 || self.EFAdj(eid, 2) < 1
                be = true;
            else
                be = false;
            end
        end
        
        function bv = isBoundaryVertex(self, vid)
            for eid = self.VEAdj
                if self.isBoundaryEdge(eid)
                    bv = true;
                    return
                end
            end
            bv = false;
        end
        
        function edge_exists = edgeExists(self, eid)
            error('Not implemented')
            if sum(self.EFAdj(eid,:)) > 0
                edge_exists = true;
            else
                edge_exists = false;
            end
        end
        
        function res = oriented_edge_is_in_face(self, fid, eid)
            v1 = self.EVAdj(eid, 1);
            v2 = self.EVAdj(eid, 2);
            verts = self.F(fid, :);
            if ( v1 == verts(1) && v2 == verts(2) ) || ...
               ( v1 == verts(2) && v2 == verts(3) ) || ...
               ( v1 == verts(3) && v2 == verts(1) )
                res = true;
            else
                res = false;
            end
        end

        function res = directed_edge_is_in_face(self, fid, vid1, vid2)
            verts = self.F(fid, :);
            if ( vid1 == verts(1) && vid2 == verts(2) ) || ...
               ( vid1 == verts(2) && vid2 == verts(3) ) || ...
               ( vid1 == verts(3) && vid2 == verts(1) )
                res = true;
            else
                res = false;
            end
        end
        
        % General info (number of edges, genus, etc)
        
        function nbE = get.nbE(self)
            if isempty(self.cached_nbE)
                %[I, ~] = find(self.EFAdj > 0);
                %I = [0; sort(I); inf];
                %self.cached_nbE = sum(I(2:end-1) > I(1:end-2) & I(2:end-1) < I(3:end));
                %self.cached_nE = sum(I(2:end-1) > I(1:end-2));
                %self.cached_niE = self.cached_nE - self.cached_nbE;
                self.compute_edge_topology();
            end
            nbE = self.cached_nbE;
        end
        
        function nE = get.nE(self)
            if isempty(self.cached_nE)
                %[I, ~] = find(self.EFAdj > 0);
                %I = [0; sort(I); inf];
                %self.cached_nbE = sum(I(2:end-1) > I(1:end-2) ...
                %    & I(2:end-1) < I(3:end));
                %self.cached_nE = sum(I(2:end-1) > I(1:end-2));
                %self.cached_niE = self.cached_nE - self.cached_nbE;
                self.compute_edge_topology();
            end
            nE = self.cached_nE;
        end
        
        function niE = get.niE(self)
            if isempty(self.cached_niE)
                %[I, ~] = find(self.EFAdj > 0);
                %I = [0; sort(I); inf];
                %self.cached_nbE = sum(I(2:end-1) > I(1:end-2) & I(2:end-1) < I(3:end));
                %self.cached_nE = sum(I(2:end-1) > I(1:end-2));
                %self.cached_niE = self.cached_nE - self.cached_nbE;
                self.compute_edge_topology();
            end
            niE = self.cached_niE;
        end

        function avg = get.avg_length(self)
            %error('broken')
            %[I, ~] = find(self.EFAdj > 0);
            %I = unique(I);
            %assert(length(I) == self.nE);
            %V1 = ceil(I / self.nV);
            %V2 = mod(I-1, self.nV) + 1;
            V1 = self.EVAdj(:, 1);
            V2 = self.EVAdj(:, 2);
            avg = mean(sqrt(sum((self.V(V1,:)-self.V(V2,:)).^2, 2)));
        end
        
        function genus = get.genus(self)
            genus = round(1 - sum(self.gaussian_cur)/(4*pi));
        end
        
        function valence = getValence(obj)
            valence = sum(obj.VVAdj, 2) + sum(obj.VVAdj, 1)';
        end
        
        % Normals
        
        function FNormals = get.FNormals(self)
            if isempty(self.cached_FNormals)
                V1 = self.V(self.F(:,2), :) - self.V(self.F(:,1), :);
                V2 = self.V(self.F(:,3), :) - self.V(self.F(:,1), :);
                self.cached_FNormals = cross(V1, V2, 2);
                self.cached_FNormals = ...
                    self.cached_FNormals ./ ...
                    repmat(row_norm(self.cached_FNormals), 1, 3);
                
                self.cached_FNormals(isnan(self.cached_FNormals)) = 0;
            end
            FNormals = self.cached_FNormals;
        end
        
        function VNormals = get.VNormals(self)
            if isempty(self.cached_VNormals)
                I1 = self.F(:, 1);
                I2 = self.F(:, 2);
                I3 = self.F(:, 3);
                J = (1:self.nF)';
                S = repmat(self.AF, 3, 1); %ones(3*self.nF, 1);
                vfadj = sparse(...
                    [I1; I2; I3], ...
                    [J; J; J], ...
                    S, ...
                    self.nV, ...
                    self.nF);
                self.cached_VNormals = vfadj * self.FNormals;
                
                % normalize
                self.cached_VNormals = ...
                    self.cached_VNormals ./ ...
                    repmat(row_norm(self.cached_VNormals), 1, 3);
            end
            VNormals = self.cached_VNormals;
        end
        
        function VF_angles = get.VF_angles(self)
            if isempty(self.cached_VF_angles)
                E12 = normalize_rows(self.V(self.F(:,2),:)-self.V(self.F(:,1),:));
                E13 = normalize_rows(self.V(self.F(:,3),:)-self.V(self.F(:,1),:));
                E21 = normalize_rows(self.V(self.F(:,1),:)-self.V(self.F(:,2),:));
                E23 = normalize_rows(self.V(self.F(:,3),:)-self.V(self.F(:,2),:));
                E31 = normalize_rows(self.V(self.F(:,1),:)-self.V(self.F(:,3),:));
                E32 = normalize_rows(self.V(self.F(:,2),:)-self.V(self.F(:,3),:));
                theta1 = acos(dot(E12, E13, 2));
                theta2 = acos(dot(E23, E21, 2));
                theta3 = acos(dot(E31, E32, 2));    

                I1 = self.F(:, 1);
                I2 = self.F(:, 2);
                I3 = self.F(:, 3);
                J = (1:self.nF)';

                self.cached_VF_angles = sparse([I1;I2;I3], ...
                                               [J;J;J], ... 
                                               [theta1;theta2;theta3], ...
                                               self.nV, ...
                                               self.nF);
            end
            VF_angles = self.cached_VF_angles;
        end
        
        % Curvature
        
        function Kg = get.gaussian_cur(self)  
            if isempty(self.cached_Kg)
                self.cached_Kg = get_gaussian_curvature(self);
            end
            Kg = self.cached_Kg;
            %Kg = get_gaussian_curvature(self);
        end
        
        function gen_def = get.generator_defects(self)
            if isempty(self.cached_generator_defects)
                self.cached_generator_defects = wrapToPi(generator_angle_defects(self));
            end
            gen_def = self.cached_generator_defects;
        end
        
        function be = get.boundaryEdge(self)
            if isempty(self.cached_boundaryEdge)
                self.cached_boundaryEdge = zeros(self.nE, 1);
                for eid = 1:self.nE
                    if self.EFAdj(eid, 1) < 1 || self.EFAdj(eid, 2) < 1
                        self.cached_boundaryEdge(eid) = true;
                    else
                        self.cached_boundaryEdge(eid) = false;
                    end
                end
            end
            be = self.cached_boundaryEdge;
        end
        
        function Mk = get.mean_cur(self)
            [L, ~] = lap(self);
            Mk = 0.5*row_norm(L*self.V);
        end
        
        function Ad = cycle_defects(self, cycles, verbose)
            if nargin < 3
                verbose = false;
            end
            F = self.F; V = self.V;
            %cycles = mesh.generator_cycles;
            n_cycles = numel(cycles);
            %assert(n_cycles == self.genus*2);
            if n_cycles == 0
                Ad = [];
                return;
            end
            Ad = zeros(n_cycles, 1);
            for i = 1:n_cycles
                cy = cycles{i}(1:end-1);
                n = length(cy);
                if verbose, disp(['loop_size = ', num2str(n)]); end
                for j = 1:n
                    fid1 = cy(j);
                    fid2 = cy(mod(j-1+1, n)+1);
                    fid3 = cy(mod(j-1+2, n)+1);
                    face1 = F(fid1, :);
                    face2 = F(fid2, :);
                    face3 = F(fid3, :);

                    common_edge = intersect(face1, face2);
                    k = 1;
                    while 1
                        v1 = face2(k);
                        v2 = face2(mod(k,3)+1);
                        if length(intersect([v1, v2], common_edge)) == 2
                            break;
                        end
                        k = k + 1;
                    end

                    v12 = [face2(k), face2(mod(k,3)+1)];
                    k = mod(k, 3) + 1;
                    v23 = [face2(k), face2(mod(k,3)+1)];
                    k = mod(k, 3) + 1;
                    v31 = [face2(k), face2(mod(k,3)+1)];

                    e1 = V(v12(2), :) - V(v12(1), :);
                    if length(intersect(v23, face3)) == 2
                        e2 = V(v23(2), :) - V(v23(1), :);
                        denom = norm(e1)*norm(e2);
                        cos_a = dot(-e1, e2) / denom;
                        if verbose, disp(['-= ', num2str(acos(cos_a))]); end
                        Ad(i) = Ad(i) - acos(cos_a);
                        %Ad(i) = Ad(i) - vec_vec_angle(-e1, e2);
                    elseif length(intersect(v31, face3)) == 2
                        e3 = V(v31(2), :) - V(v31(1), :);
                        denom = norm(e1)*norm(e3);
                        cos_a = dot(e1, -e3) / denom;
                        if verbose, disp(['+= ', num2str(acos(cos_a))]); end
                        Ad(i) = Ad(i) + acos(cos_a);
                        %Ad(i) = Ad(i) + vec_vec_angle(e1, -e3);
                    else
                        error('Bad loop')
                    end
                end
            end
            
            Ad = wrapToPi(Ad);
        end
        
        % Cycles / trees / generators
        
        function alpha = cycle_to_1form(self, cy)
            alpha = zeros(self.ne, 1);
            for j = 1:length(path) - 1
                f1 = cy(j);
                f2 = cy(j + 1);

                edge_sign = sign(self.FFo(f1, f2));
                eid = abs(self.FFo(f1, f2));
                alpha(eid) = edge_sign;
            end
        end
        
        function cycles = vert_cycles(self)
            % Construct a 1-ring cycle (face list) around each vertex.
            % Unlike vf_1ring, the faces are in a consecutive order.
            if self.nbE > 0
                error('Boundary edges are not supported')
            end
            
            vf_1ring = self.vf_1ring;
            FF = self.FFAdj;
            F  = self.F;
            cycles = cell(self.nV, 1);
            for vi = 1:self.nV
                f0 = vf_1ring{vi}(1);
                cy = [f0];
                fj = -1;
                
                % Find edge that's going into vi
                for k = 1:3
                    v1 = F(f0, k);
                    v2 = F(f0, mod(k, 3)+1);
                    if v2 == vi
                        break
                    end
                end
                vj = v1;
                
                % Walk around the 1-ring of vi
                while fj ~= f0
                    % Look at the faces adjacent to the last face in the
                    % cycle
                    faces = FF(cy(end), :);
                    for fj = faces
                        % Find the face with an edge going out of vi
                        for k = 1:3
                            v1 = F(fj, k);
                            v2 = F(fj, mod(k,3)+1);
                            if v1 == vi && v2 == vj
                                cy(end+1) = fj;
                                vj = F(fj, mod(k+1,3)+1);
                                break;
                            end
                        end
                        if cy(end) == fj
                            break
                        end
                        % If fj is in the 1-ring of vi    
                        %if ismember(vi, F(fj,:))
                        %    if length(cy) > 1 && fj == cy(end-1)
                        %        % Wrong direction, skip
                        %        continue
                        %    end
                        %    cy(end+1) = fj;
                        %    break;
                        %end
                    end
                end
                cycles{vi} = cy;
            end
        end
        
        function T = get.primalMST(self)
            T = mst(self.VVAdj);
        end
        
        function T = get.dualMST(self)
            T = mst(self.FFAdj);
        end  
        
        function cycles = get.generator_cycles(self)
            %%%%% Primal/Dual tree-co-tree decomposition. Assumes no boundaries!
            if ~isempty(self.cached_generator_cycles)
                cycles = self.cached_generator_cycles;
                return;
            end

            if nnz(self.boundaryEdge) > 0
                error('Assumes no boundaries!')
            end
            
            debug = false;

            nf = self.nF; nv = self.nV;

            % Primal graph
            s = reshape(self.F,[],1);
            t = reshape(self.F(:,[2,3,1]),[],1);
            Gp = digraph(s,t); Ap = adjacency(Gp); Gp = graph(Ap);

            ne2 = sum(sum(Ap));
            [ii,jj] = find(Ap);
            Ep = sparse(ii,jj,1:ne2,nv,nv);

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
            Xf = self.IV2F*self.V;
            %[ii,jj] = find(Adt); line([Xf(ii,1),Xf(jj,1)]',[Xf(ii,2),Xf(jj,2)]',[Xf(ii,3),Xf(jj,3)]','linewidth',2,'color','r'); 

            % num edges in dual graph not in primal spanning tree
            nednpt = sum(sum(Adnpt));
            % num edges in dual spanning tree
            nedt = sum(sum(Adt));
            % number of edges which are in neither trees = 2*g
            if ((nednpt - nedt)/2 ~= -(nv + nf - ne2/2 - 2))
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

                if debug
                    figure; self.draw('FaceAlpha', 0.9); hold on; 
                    line([Xf(f1s(i),1),Xf(f2s(i),1)]',[Xf(f1s(i),2),Xf(f2s(i),2)]',[Xf(f1s(i),3),Xf(f2s(i),3)]','linewidth',1,'color','m'); 
                end

                cycle_a = f1s(i);
                while (pred(cycle_a(end))~=0)
                    cycle_a = [cycle_a, pred(cycle_a(end))];
                    v1 = cycle_a(end-1); v2 = cycle_a(end);

                    if debug
                        line([Xf(v1,1),Xf(v2,1)]',[Xf(v1,2),Xf(v2,2)]',[Xf(v1,3),Xf(v2,3)]','linewidth',1,'color','b'); 
                    end
                end

                cycle_b = f2s(i);
                while (pred(cycle_b(end))~=0)
                    cycle_b = [cycle_b, pred(cycle_b(end))];
                    v1 = cycle_b(end-1); v2 = cycle_b(end);

                    if debug
                        line([Xf(v1,1),Xf(v2,1)]',[Xf(v1,2),Xf(v2,2)]',[Xf(v1,3),Xf(v2,3)]','linewidth',1,'color','r'); 
                    end
                end
                
                cycles{end+1} = [cycle_a(end:-1:1),cycle_b];
            end

            for i = 1:numel(cycles)
               cy = cycles{i};
               % Find duplicate nodes in cycle
               [~, unq_inds] = unique(cy);
               dup_inds = setdiff(1:numel(cy), unq_inds);
               if ~isempty(dup_inds)
                   % Use first duplicate to extract a cycle
                   % without a tail.
                   inds = find(cy == cy(dup_inds(1)));
                   cycles{i} = cy(inds(1) : inds(2));
               end
            end
            
            self.cached_generator_cycles = cycles;
        end
        
        function C = dual_cycles_to_1forms(self, cycles)
            % It's a lot faster (as of 2017a) to make a local copy than
            % call the getter each iteration. Matlab...
            FFo = self.FFo; 
            C = sparse(self.nE, numel(cycles));
            for i = 1:numel(cycles)
                cy = cycles{i};
                for j = 1:length(cy) - 1
                    f1 = cy(j);
                    f2 = cy(j + 1);

                    edge_sign = sign(FFo(f1, f2));
                    eid = abs(FFo(f1, f2));
                    C(eid, i) = edge_sign;
                end
            end
        end
        
        function C = primal_cycles_to_1forms(self, cycles)
            v1 = self.EVAdj(:, 1);
            v2 = self.EVAdj(:, 2);
            VV = sparse([v1; v2], [v2; v1], [1:self.nE, 1:self.nE]', self.nV, ...
                self.nV, 2*self.nE);
            
            C = sparse(self.nE, numel(cycles));
            for i = 1:numel(cycles)
                cy = cycles{i};
                for j = 1:length(cy) - 1
                    v1 = cy(j);
                    v2 = cy(j+1);
                    eid = VV(v1, v2);
                    C(eid, i) = sign(v2 - v1);
                end
            end
        end
        
        function H = get.H(self)
            if ~isempty(self.cached_H)
                H = self.cached_H;
                return;
            end

            cycles = self.generator_cycles();

            ng2 = numel(cycles);
            assert(ng2 == 2*self.genus)
            if ng2 < 1
                H = [];
                return
            end
            
            H = self.dual_cycles_to_1forms(cycles);

%             if debug
%                 for i = 1:ng2
%                     figure
%                     self.draw('FaceAlpha', 0.8, 'PlotField', false)
%                     %self.drawLabels()
%                     hold on
%                     cy = cycles{i};
%                     for j = 1:length(cy) - 1
%                         f1 = cy(j);
%                         f2 = cy(j+1);
%                         p1 = sum(self.V(self.F(f1, :), :), 1) ./ 3;
%                         p2 = sum(self.V(self.F(f2, :), :), 1) ./ 3;
%                         
%                         eid = abs(self.FFo(f1, f2));
%                         edge_sign = H(eid, i);
%                         
%                         if edge_sign > 0
%                             arrow(p1, p2, 'Length', 1, 'Width', 1, 'FaceColor', 'b', 'EdgeColor', 'b');
%                         else
%                             arrow(p1, p2, 'Length', 1, 'Width', 1, 'FaceColor', 'r', 'EdgeColor', 'r');
%                         end
%                     end
%                 end
%                 hold off
%             end
            
            self.cached_H = H;
        end
        
        function [Gamma, defects, cycles] = SFCBasisp(self)
            if nnz(self.boundaryEdge) > 0
                error('Assumes no boundaries!')
            end
            
            debug = false;

            nf = self.nF; nv = self.nV;

            % Primal graph
            s = reshape(self.F,[],1);
            t = reshape(self.F(:,[2,3,1]),[],1);
            Gp = digraph(s,t); Ap = adjacency(Gp); Gp = graph(Ap);

            ne2 = sum(sum(Ap));
            assert(ne2 == 2*self.nE)
            %[ii,jj] = find(Ap);
            %Ep = sparse(ii,jj,1:ne,nv,nv);
            
            % Dual spanning tree
            [Tp, pred] = minspantree(Gp,'method','sparse'); Apt = adjacency(Tp);
            
            % Dual adjacency matrix without dual spanning tree
            Apnpt = Ap - Apt;
            if (norm(Apnpt - Apnpt','fro')~=0)
                error('?');
            end
            Apnpt = Apnpt - triu(Apnpt);
            [f1s, f2s] = find(Apnpt);
            assert(length(f1s) == self.nF - 1)
            assert(self.nE - self.nF + 1 == self.nV - 1 + 2*self.genus)
            
            % Construct cycles. For each edge that's not in the dual
            % spanning tree, follow the tree back to the root from both
            % endpoints. This creates path_a and path_b, which together
            % form a cycle (with a tail).
            cycles = cell(length(f1s), 1);
            for i = 1:length(f1s)
                path_a = f1s(i);
                while (pred(path_a(end))~=0)
                    path_a = [path_a, pred(path_a(end))];
                    
                    if debug
                        v1 = path_a(end-1); v2 = path_a(end);
                        line([Xf(v1,1),Xf(v2,1)]',[Xf(v1,2),Xf(v2,2)]',[Xf(v1,3),Xf(v2,3)]','linewidth',1,'color','b'); 
                    end
                end

                path_b = f2s(i);
                while (pred(path_b(end))~=0)
                    path_b = [path_b, pred(path_b(end))];
                    
                    if debug
                        v1 = path_b(end-1); v2 = path_b(end);
                        line([Xf(v1,1),Xf(v2,1)]',[Xf(v1,2),Xf(v2,2)]',[Xf(v1,3),Xf(v2,3)]','linewidth',1,'color','r'); 
                    end
                end
                
                cycles{i} = [path_a(end:-1:1),path_b];
            end
            
            % Remove the tail from each cycle
            for i = 1:numel(cycles)
               cy = cycles{i};
               % Find duplicate nodes
               [~, unq_inds] = unique(cy);
               dup_inds = setdiff(1:numel(cy), unq_inds);
               if ~isempty(dup_inds)
                   % Take the first and second occurance of the first
                   % duplicate to form a cycle without a tail.
                   inds = find(cy == cy(dup_inds(1)));
                   cycles{i} = cy(inds(1) : inds(2));
               end
            end
            
            Gamma = self.primal_cycles_to_1forms(cycles);
            %defects = self.cycle_defects(cycles);
            defects = nan;
        end
        
        function [Gamma, defects, cycles] = SFCBasisd(self)
            % Construct a strictly fundamental cycle basis on the dual graph. 
            % Assumes no boundary.
            %
            % Output:
            %   Gamma   - ne x (nv-1+2g) matrix of basis cycles (each
            %             non-zero entry is equal to +-1 according to 
            %             the orientation)
            %   defects - the angle defect of each cycle
            %   cycles  - (nv-1+2g) x 1 cell array of cycles, each a list 
            %             of faces
            
            if nnz(self.boundaryEdge) > 0
                error('Assumes no boundaries!')
            end
            
            debug = false;

            nf = self.nF; nv = self.nV;

            % Primal graph
            s = reshape(self.F,[],1);
            t = reshape(self.F(:,[2,3,1]),[],1);
            Gp = digraph(s,t); Ap = adjacency(Gp); Gp = graph(Ap);

            ne2 = sum(sum(Ap));
            assert(ne2 == 2*self.nE)
            [ii,jj] = find(Ap);
            Ep = sparse(ii,jj,1:ne2,nv,nv);

            % Dual graph
            F1 = sparse(s,t,[1:nf,1:nf,1:nf]',nv,nv);
            F2 = sparse(t,s,[1:nf,1:nf,1:nf]',nv,nv);
            [ii,jj,ss] = find(Ep); ll = sub2ind(size(Ep),ii,jj);
            Ed = sparse(F1(ll),F2(ll),ss,nf,nf);
            Ad = double(Ed ~= 0);
            Gd = graph(Ad);
            
            % Dual spanning tree
            [Td, pred] = minspantree(Gd,'method','sparse'); Adt = adjacency(Td);
            
            % Dual adjacency matrix without dual spanning tree
            Adndt = Ad - Adt;
            if (norm(Adndt - Adndt','fro')~=0)
                error('?');
            end
            Adndt = Adndt - triu(Adndt);
            [f1s, f2s] = find(Adndt);
            assert(length(f1s) == self.nE - self.nF + 1)
            assert(self.nE - self.nF + 1 == self.nV - 1 + 2*self.genus)
            
            % Construct cycles. For each edge that's not in the dual
            % spanning tree, follow the tree back to the root from both
            % endpoints. This creates path_a and path_b, which together
            % form a cycle (with a tail).
            cycles = cell(length(f1s), 1);
            for i = 1:length(f1s)
                path_a = f1s(i);
                while (pred(path_a(end))~=0)
                    path_a = [path_a, pred(path_a(end))];
                    
                    if debug
                        v1 = path_a(end-1); v2 = path_a(end);
                        line([Xf(v1,1),Xf(v2,1)]',[Xf(v1,2),Xf(v2,2)]',[Xf(v1,3),Xf(v2,3)]','linewidth',1,'color','b'); 
                    end
                end

                path_b = f2s(i);
                while (pred(path_b(end))~=0)
                    path_b = [path_b, pred(path_b(end))];
                    
                    if debug
                        v1 = path_b(end-1); v2 = path_b(end);
                        line([Xf(v1,1),Xf(v2,1)]',[Xf(v1,2),Xf(v2,2)]',[Xf(v1,3),Xf(v2,3)]','linewidth',1,'color','r'); 
                    end
                end
                
                cycles{i} = [path_a(end:-1:1),path_b];
            end
            
            % Remove the tail from each cycle
            for i = 1:numel(cycles)
               cy = cycles{i};
               % Find duplicate nodes
               [~, unq_inds] = unique(cy);
               dup_inds = setdiff(1:numel(cy), unq_inds);
               if ~isempty(dup_inds)
                   % Take the first and second occurance of the first
                   % duplicate to form a cycle without a tail.
                   inds = find(cy == cy(dup_inds(1)));
                   cycles{i} = cy(inds(1) : inds(2));
               end
            end
            
            Gamma = self.dual_cycles_to_1forms(cycles);
            defects = self.cycle_defects(cycles);
        end
        
        function set_ffield(self, ...
                degree, ...
                local_frames, ...
                frame_diffs, ...
                ffield_angles, ...
                ffield_vectors, ...
                vert_sing, ...
                gen_sing, ...
                miq_energy, ...
                connection)
            if nargin < 10
                connection = [];
            end
            
            self.degree          = degree;
            self.local_frames    = local_frames;
            self.frame_diffs     = frame_diffs;
            self.ffield_angles   = ffield_angles;
            self.ffield_vectors  = ffield_vectors;
            self.n_vert_sing     = size(vert_sing, 1);
            self.n_gen_sing      = size(gen_sing, 1);
            self.n_singularities = size(vert_sing, 1); %+ size(gen_sing, 1);
            self.vert_sing       = vert_sing;
            self.gen_sing        = gen_sing;
            self.miq_energy      = miq_energy;
            self.connection      = connection;
        end
        
        % Other
        
        function vert_vects = face_vect_to_vertex_vect(self, face_vects)
            % Takes a vector field defined on the faces and returns
            % a vector field defined on the vertices. Each such vector
            % is computed as the weighted angle average.
            [r, ~] = size(face_vects);
            if r ~= self.nF
                face_vects = reshape(face_vects, self.nF, []);
                assert(size(face_vects, 2) == 3)
            end
                
            vert_vects = full(self.VF_angles * face_vects);
            vert_vects = vert_vects ./ repmat(full(sum(self.VF_angles, 2)), 1, 3);
        end
        
        %%%%%%%%%%%%
        % Plotting %
        %%%%%%%%%%%%

        function colorFaces(self, faces, color)
            if nargin < 3
                color = 'r';
            end
            patch('faces', self.F(faces, :), ...
                  'vertices', self.V, ...
                  'FaceColor', color, ...
                  'FaceAlpha', 0.8, ...
                  'EdgeAlpha', 0.2);
        end

        function colorEdges(self, edges, color)
            if nargin < 3
                color = 'm';
            end
            EV = self.EVAdj(edges, :);
            V1 = EV(:, 1);
            V2 = EV(:, 2);
            X = self.V(V1, 1);
            Y = self.V(V2, 2);
            Z = self.V(V2, 3);
            line(X', Y', Z', 'LineWidth', 2, 'color', color)
        end
        
        function A = montage(self, filename, funcs, rows, titles, varargin)
            if nargin < 4
                rows = 1;
            end
            if nargin < 5
                titles = {};
            end
            sz = size(funcs);
            
            %try
            for i = 1:sz(2)
                fh = figure;
                self.draw(funcs(:, i), varargin{:});
                set(gcf, 'WindowStyle', 'docked')
                set(gcf,'color','w');
                
                pngname = sprintf('%s_%04d.png', filename, i);
                pdfname = sprintf('%s_%04d.pdf', filename, i);
                figname = sprintf('%s_%04d.fig', filename, i);
                export_fig(pngname);
                export_fig(pdfname);
                saveas(gcf, figname, 'fig');
            end
            
            cols = ceil(sz(2) / rows);
            counter = 1;
            A = [];
            for i = 1:rows
                b = [];
                for j = 1:cols
                    pngname = sprintf('%s_%04d.png',filename,cols*(i-1)+j);
                    a = imread(pngname);
                    b = cat(2,b,a);
                    counter = counter + 1;
                    if counter == sz(2)+1
                        break
                    end
                end
                A = cat(1,A,b);
                imwrite(A, [filename, '.png']);
                %delete( sprintf('%s_*.png',filename) );
            end
            %catch ME
            %    delete( sprintf('%s_*.png',filename) );
            %    disp(getReport(ME, 'extended'));
            %end
            
            %for i = 1:sz(2)
            %    b = [];
            %    pngname = sprintf('%s_%04d.png', filename, i);
            %    if r <= cols
            %        a = imread(pngname);
            %        b = cat(2, b, a);
            %    else
            %        A = cat(1, A, b);
            %    end
            %end
            
        end
        
        function draw(obj, varargin)
            p = inputParser;
            
            addOptional(p, 'Func', ones(obj.nF, 1));           
            addOptional(p, 'FaceColor', []);
            addOptional(p, 'EdgeColor', 'none');
            addOptional(p, 'FaceAlpha', 1);
            addOptional(p, 'EdgeAlpha', 0.2);
            addOptional(p, 'scale', 1);
            % TODO assumes cross field
            addOptional(p, 'nrosy_colors', {'b', 'b', 'b', 'b'});
            addOptional(p, 'PlotField',false);
            addOptional(p, 'PlotSing', true);
            addOptional(p, 'Caxis', 'auto');
            addOptional(p, 'View', [0 1 0]);
            addOptional(p, 'Dock', 0);
            addOptional(p, 'Colorbar', false);
            addOptional(p, 'MarkerSize', 80);
            addOptional(p, 'Colormap', []);
            addOptional(p, 'Camera', []);
            
            parse(p, varargin{:});
            
            Func = p.Results.Func;
            FaceColor = p.Results.FaceColor;
            if isempty(FaceColor)
                if length(Func) == obj.nF
                    FaceColor = 'flat';
                elseif length(Func) == obj.nV
                    FaceColor = 'interp';
                else
                    warning('Given function is not of size nv or nf');
                    FaceColor = 'w';
                end
            end
            
            if size(Func, 2) == 1
                cw = 'CData';
            elseif size(f, 2) == 3
                cw = 'FaceVertexCData';
            end
            
            hold on
            
            % Plot singularities
            if ~isempty(obj.vert_sing) && p.Results.PlotSing
                S = obj.vert_sing;
                colors = linspecer(2);
                for j = 1:size(S,1)
                    fid = S(j, 1);
                    if S(j, 2) > 0
                        H = plot3( obj.V(fid,1), obj.V(fid,2), obj.V(fid,3), '.', 'color', colors(1,:) );
                    else
                        H = plot3( obj.V(fid,1), obj.V(fid,2), obj.V(fid,3), '.', 'color', colors(2,:) );
                    end
                    set( H, 'MarkerSize', p.Results.MarkerSize );
                end
            end

            % Plot nrosy field
            if ~isempty(obj.ffield_vectors) && p.Results.PlotField
                ffield = obj.ffield_vectors;
                nrosy_colors = p.Results.nrosy_colors;
                scale = p.Results.scale;
                for j = 0:obj.degree-1
                    inds = 1+j*obj.nF : obj.nF+j*obj.nF;
                    obj.drawFaceField(ffield(inds, :), ...
                        'color', nrosy_colors{j+1}, ...
                        'AutoScale', 'on', ...
                        'AutoScaleFactor', scale)
                end
            end

            patch('faces', obj.F, ...
                  'vertices', obj.V, ...
                  cw, Func, ...
                  'FaceColor', FaceColor, ...
                  'EdgeColor', p.Results.EdgeColor, ...
                  'FaceAlpha', p.Results.FaceAlpha, ...
                  'EdgeAlpha', p.Results.EdgeAlpha);
            if p.Results.Colorbar
                colorbar
            end
                  
              
            hold off

            ax = gca;
            ax.Clipping = 'off';
            
            cameratoolbar; cameratoolbar('SetCoordSys','none'); 
            axis equal; axis off;
            
            if ~isempty(p.Results.View)
                view(p.Results.View);
            end
            if ~isempty(p.Results.Camera)
                set_camera(gca, p.Results.Camera);
            end
            
            camlight('headlight')
            lighting phong
            material dull
            
            caxis(p.Results.Caxis);
            if ~isempty(p.Results.Colormap)
                colormap(p.Results.Colormap);
            end
            if p.Results.Dock
                set(gcf, 'WindowStyle', 'docked');
            end
        end
        
        function drawFaceField(self, vectors, varargin)
            P1 = (self.V(self.F(:,1), :) + ...
                  self.V(self.F(:,2), :) + ...
                  self.V(self.F(:,3), :)) / 3;
            px = P1(:,1);
            py = P1(:,2);
            pz = P1(:,3);
            u = vectors(:,1);
            v = vectors(:,2);
            w = vectors(:,3);
            quiver3(px, py, pz, u, v, w, varargin{:});
        end
        
        function drawField(self, vectors, varargin)
            [r, c] = size(vectors);
            if c == 1
                vectors = reshape(vectors, [], 3);
                [r, c] = size(vectors);
            end
            %vectors(row_norm(vectors)<1e-3) = 0;
            if r == self.nV
                px = self.V(:,1);
                py = self.V(:,2);
                pz = self.V(:,3);
            elseif r == self.nF
                P1 = (self.V(self.F(:,1), :) ...
                      + self.V(self.F(:,2), :) ...
                      + self.V(self.F(:,3), :)) / 3;
                px = P1(:,1);
                py = P1(:,2);
                pz = P1(:,3);
            end
            %vectors = vectors ./ min(row_norm(vectors));
            %vectors = normalize_rows(vectors);
            %vectors = vectors / max(row_norm(vectors));
            u = vectors(:,1);
            v = vectors(:,2);
            w = vectors(:,3);
            quiver3(px, py, pz, u, v, w, varargin{:});
        end
        
        function drawFNormals(self, varargin)
            hold on
            P1 = (self.V(self.F(:,1), :) ...
                + self.V(self.F(:,2), :) ...
                + self.V(self.F(:,3), :)) / 3;
            x = P1(:,1);
            y = P1(:,2);
            z = P1(:,3);
            u = self.FNormals(:,1);
            v = self.FNormals(:,2);
            w = self.FNormals(:,3);
            quiver3(x, y, z, u, v, w, varargin{:})
        end
        
        function drawVNormals(self, varargin)
            x = self.V(:,1);
            y = self.V(:,2);
            z = self.V(:,3);
            u = self.VNormals(:,1);
            v = self.VNormals(:,2);
            w = self.VNormals(:,3);
            quiver3(x, y, z, u, v, w, varargin{:});
        end
        
        function drawLabels(self, varargin)
            p = inputParser;
            %addOptional(p, 'function', zeros(m.nF, 1));
            
            addOptional(p, 'Vertex', true);           
            addOptional(p, 'Face', true);
            addOptional(p, 'Edge', true);
            addOptional(p, 'FontSize', 14);

            parse(p, varargin{:});
            vertex = p.Results.Vertex;
            face = p.Results.Face;
            edge = p.Results.Edge;
            fnt_sz = p.Results.FontSize;
            
            if vertex
                for vid = 1:self.nV
                    p = self.V(vid,:);
                    p = p + 0.001*self.VNormals(vid, :);
                    text(p(1), p(2), p(3), ['v', num2str(vid)], 'FontSize', fnt_sz)
                end
            end

            if face
                for fid = 1:self.nF
                    vid1 = self.F(fid, 1);
                    vid2 = self.F(fid, 2);
                    vid3 = self.F(fid, 3);
                    p = (self.V(vid1,:)+self.V(vid2,:)+self.V(vid3,:)) / 3;
                    p = p + 0.1*self.FNormals(fid, :);
                    text(p(1), p(2), p(3), ['f', num2str(fid)], 'FontSize', fnt_sz)
                end
            end
            
            if edge
                EF = self.EFAdj;
                for eid = 1:self.nE
                    vid1 = self.EVAdj(eid, 1);
                    vid2 = self.EVAdj(eid, 2);
                    p = (self.V(vid1,:) + self.V(vid2,:)) ./ 2;
                    text(p(1), p(2), p(3), ['e', num2str(eid)], 'FontSize', fnt_sz)
                    
                    p = p + 0.1 * (self.V(vid2,:) - self.V(vid1,:));
                    f1 = EF(eid, 1);
                    f2 = EF(eid, 2);
                    if f1 >0 && f2 > 0
                        fn1 = self.FNormals(f1, :);
                        fn2 = self.FNormals(f2, :);
                        p = p + 0.01 * (fn1 + fn2) / 2; 
                        text(p(1), p(2), p(3), ['e', num2str(eid)], 'FontSize', fnt_sz)
                    end
                end
            end
        end
        
        function labelFaces(self, f_labels, FORMAT, offset, varargin)
            if nargin < 2
                f_labels = 1:self.nF;
            end
            if nargin < 3
                FORMAT = '%d';
            end
            if nargin < 4
                offset = [0, 0, 0];
            end
            for fid = 1:self.nF
                v1 = self.V(self.F(fid, 1), :);
                v2 = self.V(self.F(fid, 2), :);
                v3 = self.V(self.F(fid, 3), :);
                p = (v1 + v2 + v3) / 3 + offset;
                text(p(1), p(2), p(3), num2str(f_labels(fid), FORMAT), varargin{:});
            end
        end
        
        function labelEdges(self, inds, labels, FORMAT, offset, varargin)
            if nargin < 2
                inds = 1:self.nE;
            end
            if nargin < 3
                labels = 1:self.nE;
            end
            if nargin < 4
                FORMAT = '%d';
            end
            if nargin < 5
                offset = [0, 0, 0];
            end
            for i = 1:length(inds)
                eid = inds(i);
                
                v1 = self.V(self.EVAdj(eid, 1), :);
                v2 = self.V(self.EVAdj(eid, 2), :);
                p = (v1 + v2) / 2 + offset;
                
                fid1 = self.EFAdj(eid, 1);
                fid2 = self.EFAdj(eid, 2);
                if fid1 == 0
                    fid1 = fid2;
                end
                if fid2 == 0
                    fid2 = fid1;
                end
                normal = (self.FNormals(fid1, :) + self.FNormals(fid2, :)) / 2;
                p = p + 0.1 * self.avg_length * normal;
                
                text(p(1), p(2), p(3), num2str(labels(i), FORMAT), varargin{:});
            end
        end
        
        function labelVertices(self, v_labels, FORMAT, offset, varargin)
            if nargin < 2
                v_labels = 1:self.nV;
            end
            if nargin < 3
                FORMAT = '%d';
            end
            if nargin < 4
                offset = [0, 0, 0];
            end
            for vid = 1:self.nV
                p = self.V(vid, :) + offset;
                text(p(1), p(2), p(3), num2str(v_labels(vid), FORMAT), varargin{:});
            end
        end
        
        function plotCycles(self, cycles)
            col = linspecer(numel(cycles));
            self.draw('FaceColor', 'w', 'FaceAlpha', 1); hold on
            EF = self.EFAdj;
            FFo = self.FFo;
            % The midpoint of each face, translated a bit in the normal
            % direction.
            Xf = self.IV2F * self.V + 0.1*self.avg_length*self.FNormals;
            axis(axis)
            for i = 1:numel(cycles)
                cy = cycles{i};
                t = 0.1*self.avg_length*randn;
                Xf = Xf + t * self.FNormals;
                for j = 1:length(cy)-1
                    f1 = cy(j);
                    f2 = cy(j+1);
                    eid = abs(FFo(f1, f2));
                    %fn1 = self.FNormals(f1, :);
                    %fn2 = self.FNormals(f2, :);
                    %p1 = Xf(f1, :);
                    %p2 = Xf(f2, :);
                    %p1 = sum(self.V(self.F(f1, :), :), 1) ./ 3;
                    %p2 = sum(self.V(self.F(f2, :), :), 1) ./ 3;
                    
                    %arrow(p1, p2, 'Length', 5, 'EdgeColor', col(i, :));
                    line([Xf(f1,1),Xf(f2,1)]', ...
                        [Xf(f1,2),Xf(f2,2)]', ...
                        [Xf(f1,3),Xf(f2,3)]', ...
                        'linewidth',1,'color',col(i,:)); 
                    
                    %edge_sign = Gamma(eid, i);
%                     p = (0.7*p1 + 0.3*p2);
%                     if edge_sign > 0
%                         text(p(1), p(2), p(3), '+', 'color', 'b', 'FontSize', 14);
%                     else
%                         text(p(1), p(2), p(3), '-', 'color', 'r', 'FontSize', 14);
%                     end
                end
                Xf = Xf - t * self.FNormals;
            end
            hold off
        end
        
        function plotH(self, inds)
            colors = {'r', 'b', 'm', 'c', 'g', 'k'};
            
            ng = size(self.H, 2);
            if nargin < 2
                inds = 1:ng;
            end
            
            EF = self.EFAdj;
            
            for i = inds
                %figure(); hold on; self.draw('FaceAlpha', 0.8);
                %fc = 1 / i;
                %colorID = max(1, sum(fc > [0:1/length(cm(:,1)):1])); 
                %myColor = cm(colorID, :); % returns your color
                if i <= length(colors)
                    myColor = colors{i};
                else
                    myColor = 'r';
                end
                
                eids = find(self.H(:, i) ~= 0);
                for j = 1:length(eids)
                    eid = eids(j);
                    f1 = EF(eid, 1);
                    f2 = EF(eid, 2);
                    p1 = sum(self.V(self.F(f1, :), :), 1) ./ 3;
                    p2 = sum(self.V(self.F(f2, :), :), 1) ./ 3;
                    edge_sign = self.H(eid, i);
                    
                    arrow(p1, p2, 'Length', 5, 'EdgeColor', myColor);
                    p = (0.7*p1 + 0.3*p2);
                    if edge_sign > 0
                        text(p(1), p(2), p(3), '+', 'color', 'b', 'FontSize', 14);
                    else
                        text(p(1), p(2), p(3), '-', 'color', 'r', 'FontSize', 14);
                    end
                end
            end
            %hold off
        end
        
        function plotEdgePaths(self, paths)
            colors = {'r', 'b', 'm', 'c', 'g', 'k'};
            EF = self.EFAdj;
            
            for i = 1:numel(paths)
                %path = paths{i};
                if i <= length(colors)
                    myColor = colors{i};
                else
                    myColor = 'r';
                end
                
                %eids = find(self.H(:, i) ~= 0);
                path = paths{i};
                for j = 1:length(path)
                    eid = path(j);
                    f1 = EF(eid, 1);
                    f2 = EF(eid, 2);
                    p1 = sum(self.V(self.F(f1, :), :), 1) ./ 3;
                    p2 = sum(self.V(self.F(f2, :), :), 1) ./ 3;
                    %edge_sign = self.H(eid, i);
                    
                    arrow(p1, p2, 'Length', 1, 'EdgeColor', myColor);
                    p = (0.7*p1 + 0.3*p2);
                    if edge_sign > 0
                        text(p(1), p(2), p(3), 'pos', 'color', 'b');
                    else
                        text(p(1), p(2), p(3), 'neg', 'color', 'r');
                    end
                end
            end
        end
        
    end
end

