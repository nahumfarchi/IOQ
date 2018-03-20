function [period_is_fixed, period] = reduce_search_space(mesh, face_is_constrained, thetas, k)
    %REDUCE_SEARCH_SPACE

    import java.util.LinkedList

    EF = mesh.EFAdj;
    FE = mesh.FEAdj;
    FF = mesh.FFAdj;
    period_is_fixed = zeros(mesh.nE, 1);
    period = zeros(mesh.nE, 1);

    start = zeros(mesh.nF, 1);
    q = LinkedList();
    for fid = 1:mesh.nF
        if face_is_constrained(fid)
            %q.add(fid);
            q.addLast(fid);
            start(fid) = true;
        end
    end

    visited = zeros(mesh.nF, 1);

    while ~q.isEmpty()
        %current_fid = q.remove();
        current_fid = q.removeFirst();
        visited(current_fid) = true;

        for i = 1:3
            neighbor_eid = FE(current_fid, i);
            neighbor_fid = FF(current_fid, i);
            if neighbor_fid >= 1
                if ~visited(neighbor_fid) && ~start(neighbor_fid)
                    period_is_fixed(neighbor_eid) = true;
                    period(neighbor_eid) = 0;
                    visited(neighbor_fid) = true;
                    %q.push(neighbor_fid);
                    q.addLast(neighbor_fid)
                end
            else
                period_is_fixed(neighbor_eid) = true;
                period(neighbor_eid) = 0;
            end
        end
    end

    for eid = 1:mesh.nE
        fid1 = EF(eid, 1);
        fid2 = EF(eid, 2);

        if fid1 >= 1 && fid2 >= 1 && face_is_constrained(fid1) && face_is_constrained(fid2)
            period_is_fixed(eid) = true;
            period(eid) = 2.0 / pi * (thetas(fid2) - thetas(fid1) - k(eid));
        end
    end

end

