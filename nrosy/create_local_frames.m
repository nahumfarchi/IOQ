function [local_frames, r] = create_local_frames(mesh)
    % function [local_frames, r] = create_local_frames(mesh)
    %
    % Create orthogonal local frames for each face in the given mesh.
    %
    % Output:
    %    local_frames - 2*nF x 3, first nF rows are for e1, last nF rows are for
    %                   e2.
    %    r - nEx1 vector of frame differences according to orientation, that is, angle(e_right) - angle(e_left).

    e1 = mesh.V(mesh.F(:,2),:) - mesh.V(mesh.F(:,1),:);
    e1 = normalize_rows(e1);
    e2 = normalize_rows(cross(mesh.FNormals, e1, 2));
    local_frames = [e1; e2];
    
    % Calculate frame differences
    if nargout > 1
        r = calc_diffs(mesh, local_frames);
    end

end

function r = calc_diffs(mesh, local_frames)
    V = mesh.V;
    F = mesh.F;
    EF = mesh.EFAdj;
    EV = mesh.EVAdj;
    nE = mesh.nE;
    
    r = zeros(nE, 1);
    
    be = mesh.boundaryEdge;
    for eid = 1:nE
        if be(eid)
            continue
        end
        % triangle abc
        fi = EF(eid, 1);
        % triangle abd
        fj = EF(eid, 2);
        
        a = V(EV(eid, 1), :)';
        b = V(EV(eid, 2), :)';
        if F(fi,1) ~= EV(eid,1) && F(fi,1) ~= EV(eid,2)
            cid = F(fi,1);
        elseif F(fi,2) ~= EV(eid,1) && F(fi,2) ~= EV(eid,2)
            cid = F(fi,2);
        elseif F(fi,3) ~= EV(eid,1) && F(fi,3) ~= EV(eid,2)
            cid = F(fi,3);
        end
        if F(fj,1) ~= EV(eid,1) && F(fj,1) ~= EV(eid,2)
            did = F(fj,1);
        elseif F(fj,2) ~= EV(eid,1) && F(fj,2) ~= EV(eid,2)
            did = F(fj,2);
        elseif F(fj,3) ~= EV(eid,1) && F(fj,3) ~= EV(eid,2)
            did = F(fj,3);
        end
        
        c = V(cid, :)';
        d = V(did, :)';
        %c = V(setdiff(F(fi, :), EV(eid,:)), :)';
        %d = V(setdiff(F(fj, :), EV(eid, :)), :)';
        
        %   a
        % c   d
        %   b
        ab = b - a; % common edge
        ac = c - a;
        ad = d - a;
        
        Ei = ortho(ab, ac);
        Ej = ortho(ab, -ad);
        
        % Transfer each frame vector to its local basis
        e1i = Ei' * local_frames(fi, :)';
        e1j = Ej' * local_frames(fj, :)';
        
        % Calculate the angle between them
        r(eid) = -(atan2(e1j(2), e1j(1)) - atan2(e1i(2), e1i(1)));
    end
end

function E = ortho(u, v)
    u = u / norm(u);
    v = v - (v'*u) * u;
    v = v / norm(v);
    w = cross(u, v);
    %w = w / norm(w);
    E = [u, v, w];
end

