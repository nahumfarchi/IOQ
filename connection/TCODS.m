function [resm] = TCODS(m, varargin)
    % function [m] = TCODS(m, varargin)
    %
    % Create a nrosy field by finding the trivial connection that is closest
    % to the Levi-Civita connection.
    %
    % Input:
    %   TODO outdated
    %   m - mesh object or file path (if using mex this has to be a path)
    %   S - nSx2 matrix of (fid,ki) pairs, where ki is the singularity
    %       index
    %   degree - field degree
    %   f0 - starting face
    %   theta0 - starting angle
    %   verbose - true/false
    %
    % Output:
    %   res - a copy of the input mesh m with the field
    %
    % Example:
    %   TCODS(m, ...
    %         'k', [alpha;beta], ...
    %         'f0', THETA0, ...
    %         'theta0', THETA0, ...
    %         'degree', DEGREE, ...
    %         'CreateFField', true, ...
    %         'Duplicate', true, ...
    %         'gConstraintVec', GVEC);
    
    p = inputParser;
    addOptional(p, 'sing', []);
    addOptional(p, 'k', []);
    addOptional(p, 'f0', 1);
    addOptional(p, 'theta0', 0);
    addOptional(p, 'BoundaryConstraints', []);
    addOptional(p, 'gConstraintVec', []);
    addOptional(p, 'degree', 4);
    addOptional(p, 'verbose', true);
    addOptional(p, 'log', -1);
    addOptional(p, 'CreateFField', true);
    addOptional(p, 'Duplicate', true);
    addOptional(p, 'Constraints', []);
    addOptional(p, 'Gamma', []);
    
    parse(p, varargin{:});
    opt = p.Results;
    
    V = m.V; F = m.F;
    nv = m.nV; nf = m.nF; ng2 = 2*m.genus;
    sing = opt.sing;
    k = opt.k;
    if isempty(k) && isempty(sing)
        error('TCODS requires atleast one of the arguments k or sing')
    elseif ~isempty(sing) && isempty(k)
        k = zeros(nv+ng2, 1);
        k(sing(:,1)) = sing(:,2);
    else
        sing = find(k); sing = [sing, k(sing)];
    end
    vert_sing = find(sing(:,1)<=nv & sing(:,2)~=0); vert_sing = [sing(vert_sing, 1), k(sing(vert_sing),1)];
    gen_sing = find(sing(:,1)>nv & sing(:,2)~=0); gen_sing = [sing(gen_sing, 1), k(sing(gen_sing, 1))];
    f0 = opt.f0;
    theta0 = opt.theta0;
    gvec = opt.gConstraintVec;
    degree = opt.degree;
    verbose = opt.verbose;
    log = opt.log;
    constraints = opt.Constraints;
    nc = size(constraints, 1);
    gamma = opt.Gamma;
    if isempty(gamma)
        gamma = zeros(max(0, nc-1), 1);
    end
    
    [local_frames, frame_diffs] = create_local_frames(m);
    
    [A, K, ~, d1, ~] = tcods_gsystem(V, F);
    %H = -m.H;
    %K = [get_gaussian_curvature(m); mod(generator_angle_defects(m),2*pi)-pi];
    %K(end-2) = K(end-2)+2*pi;
    %[d0, d1] = get_exterior_derivatives(m);
    %A = [d0'; H'];
    
    n = size(A, 2); 
    b = -K + 2 * pi * k / degree;
    
    if ~isempty(constraints)
        constrained_faces = constraints(:, 1);
        constraint_thetas = constraints(:, 2);
        f0 = constraints(1, 1); 
        theta0 = constraints(1, 2);
        frame1 = local_frames(f0, :); 
        frame2 = local_frames(f0+nf);
        gvec = cos(theta0)*frame1 + sin(theta0)*frame2;
        [C, bc] = create_constraints_mat(m, constrained_faces, constraint_thetas, frame_diffs);
        A = [A; C];
        b = [b; pi/2 * gamma - bc];
    end
    
    x = lsqlin(speye(n, n), zeros(n, 1), [], [], A, b, -inf(n, 1), inf(n, 1));
    E = norm(x)^2;
    
    if verbose
        check_norm('A*x', 'b', 'Log', log);
        check_norm('norm(d1*x)', '0', 'Log', log);
    end
    
    if opt.CreateFField
        [ffield, theta, local_frames, frame_diffs] = connection_to_nrosy(...
                     m, ...
                     x, ...
                     f0(1), ...
                    theta0(1), ...
                    degree, ...
                    'gConstraintVec', gvec, ...
                    'LocalFrames', local_frames, ...
                    'FrameDiffs', frame_diffs);
    else
        ffield = [];
        theta = [];
        local_frames = [];
        frame_diffs = [];
    end
    
    if opt.Duplicate
        resm = copy(m);
    else
        resm = m;
    end
    resm.set_ffield(degree, ...
                    local_frames, ...
                    frame_diffs, ...
                    theta, ...
                    ffield, ...
                    vert_sing, ...
                    gen_sing, ...
                    E, ...
                    x);
    
    %S(:, 2) = S(:, 2) / degree;
    
%     p = inputParser;
%     addOptional(p, 'Verbose', true);
%     addOptional(p, 'BoundaryConstraints', []);
% 	addOptional(p, 'Mex', false);
%     % Constraint on f0 specified in global coordinates
%     addOptional(p, 'gConstraintVec', []);
%     parse(p, varargin{:});
%     
% 	if p.Results.Mex
%         if ~isempty(p.Results.gConstraintVec)
%             error('Not implemented!')
%         end
%         if isempty(gen_sing)
%             [ffield, E, connection] = ...
%                 tcods_mex(m, vert_sing(:,1), vert_sing(:,2), [], [], f0, theta0, degree);
%         else
%             [ffield, E, connection] = ...
%                 tcods_mex(m, vert_sing(:,1), vert_sing(:,2), gen_sing(:,1), gen_sing(:,2), f0, theta0, degree);
%         end
%         
%         check_norm('E', 'norm(connection)^2', 'Log', -1);
%             
%         m = Mesh(m);
%         theta = nan(m.nF, 1); % TODO
%         [local_frames, frame_diffs] = create_local_frames(m);
% 	else
% 		[connection, E] = create_trivial_connection(m, ...
% 			vert_sing, ...
%             gen_sing, ...
% 			f0, ...
% 			theta0, ...
% 			p.Results.BoundaryConstraints, ...
% 			p.Results.Verbose);
%         [ffield, theta, local_frames, frame_diffs] = connection_to_nrosy(...
%                 m, ...
%                 connection, ...
%                 f0(1), ...
%                 theta0(1), ...
%                 degree, ...
%                 'gConstraintVec', p.Results.gConstraintVec);
%         %nrosy.ffield = ffield;
%         %nrosy.theta = theta;
%         %nrosy.local_frames = local_frames;
%         %nrosy.frame_diffs = frame_diffs;
% 	end
%     
%     %nrosy.degree = degree;
%     %nrosy.S = S;
%     
%     m.set_ffield(degree, ...
%                  local_frames, ...
%                  frame_diffs, ...
%                  theta, ...
%                  ffield, ...
%                  vert_sing, ...
%                  gen_sing, ...
%                  E, ...
%                  connection);
    
end

