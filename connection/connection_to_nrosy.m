function [ffield, thetas, local_frames, frame_diffs] = connection_to_nrosy(...
    m, ...
    connection, ...
    f0, ...
    theta0, ...
    N, ...
    varargin)
    % function [ffield, thetas] = connection_to_nrosy(m, connection, f0, theta0, N, varargin)
    % Create a direction field on the mesh from the given connection, starting
    % from f0 with angle theta0.
    %
    % Input:
    %   mesh object
    %   connection - E x 1 connection vector
    %   f0 - starting face
    %   theta0 - starting angle
    %   'LocalFrames', local_frames
    %   'FrameDiffs', frame_diffs
    %   'gConstraintVec', gVec
    %
    % Output:
    %   ffield
    %   thetas - F x 1 vector of angles
    %   local_frames
    %   frame_diffs
    
    p = inputParser;
    addOptional(p, 'LocalFrames', []);
    addOptional(p, 'FrameDiffs', []);
    % Constraint on f0 specified in global coordinates
    addOptional(p, 'gConstraintVec', []);
    parse(p, varargin{:});
    
    debug = false;
    
    local_frames = p.Results.LocalFrames;
    frame_diffs = p.Results.FrameDiffs;
    if isempty(local_frames ) || isempty(frame_diffs)
        [local_frames, frame_diffs] = create_local_frames(m);
    end
    
    g_constraint_vec = p.Results.gConstraintVec;
    if ~isempty(g_constraint_vec)
        frame = [local_frames(f0, :); ...
                 local_frames(f0+m.nF,:)];
        local_v = frame * g_constraint_vec(:);
        %local_v = frame' \ g_constraint_vec(:);
        theta0 = atan2(local_v(2), local_v(1));
        fprintf('theta0=%.2f\n', theta0);
    end
    
    import java.util.LinkedList
    
    EF = m.EFAdj;
    %FE = m.FEAdj;
    FF = m.FFAdj;
    
    %if nargin < 6
    %    [local_frames, frame_diffs] = create_local_frames(m);
    %end
    %frame_diffs = compute_frame_diffs(m); % TODO stop using this
    %local_frames = create_local_frames(m);   
    
    %theta0 = constraints_to_local_angles(m, local_frames, f0, v0);

    FF_to_frame_diff = sparse(m.nF, m.nF);
    FF_to_connection = sparse(m.nF, m.nF);
    be = m.boundaryEdge;
    for eid = 1:m.nE
        if be(eid)
            continue;
        end
        fi = EF(eid, 1);
        fj = EF(eid, 2);
        FF_to_frame_diff(fi, fj) = frame_diffs(eid);
        FF_to_frame_diff(fj, fi) = -frame_diffs(eid);
        FF_to_connection(fi, fj) = connection(eid);
        FF_to_connection(fj, fi) = -connection(eid);
        
        %FF_to_frame_diff(fi, fj) = -frame_diffs(eid);
        %FF_to_frame_diff(fj, fi) = frame_diffs(eid);
        %FF_to_connection(fi, fj) = -connection(eid);
        %FF_to_connection(fj, fi) = connection(eid);
        
        %FF_to_frame_diff(fi, fj) = frame_diffs(eid); +
        %FF_to_frame_diff(fj, fi) = -frame_diffs(eid); -
        %FF_to_connection(fi, fj) = -connection(eid); -
        %FF_to_connection(fj, fi) = connection(eid); +
        
        %FF_to_frame_diff(fi, fj) = -frame_diffs(eid);
        %FF_to_frame_diff(fj, fi) = frame_diffs(eid);
        %FF_to_connection(fi, fj) = connection(eid);
        %FF_to_connection(fj, fi) = -connection(eid);
    end
    
    queue = LinkedList();
    queue.addLast(f0);
    visited = zeros(m.nF, 1);
    visisted(f0) = true;
    thetas = zeros(m.nF, 1);
    thetas(f0) = theta0;
    
    if debug
        figure()
        m.draw('FaceAlpha', 0.8)
        hold on
    end

    while ~queue.isEmpty()
        fi = queue.removeFirst();

        for fj = FF(fi, :)
            if fj > 0 && ~visited(fj)
                thetas(fj) = thetas(fi) + FF_to_frame_diff(fi, fj) + FF_to_connection(fi, fj);
                visited(fj) = true;
                queue.addLast(fj);
                
                
                if debug
                    pi = sum(m.V(m.F(fi, :), :), 1) ./ 3;
                    pj = sum(m.V(m.F(fj, :), :), 1) ./ 3;
                    arrow(pi, pj, 'Length', 1, 'color', 'r')
                end
            end
        end
    end
    
    
    %DF = local_angles_to_gvectors(m, theta, local_frames, N);
    ffield = angles_to_ffield(thetas, local_frames, N);
    %theta = thetas;
    if debug
        m.ffield_vectors = ffield;
        m.degree = N;
        m.draw('FaceAlpha', 0.8)
        hold off
    end
end

