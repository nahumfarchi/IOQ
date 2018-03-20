function [C, b] = create_constraints_mat(m, fids_list, thetas_list, frame_diffs)
    % function [C, b] = create_constraints_mat(m, fids_list, thetas_list, frame_diffs)
    %
    % Create a constraint matrix by building a
    % spanning tree (on the faces) rooted at ids_list(1). Each row in the constraint
    % matrix is a path from one of ids_list(2:end) up the tree until it
    % meets some constrained face (might be the root).
    %
    % Input:
    %   m - mesh
    %   fids_list - list of constrained face ids
    %   thetas_list - corresponding constraint angles
    %
    % Output:
    %   C - constraint paths
    %   b - The RHS. Holds thetaj - thetai_PT.
    %
    % When creating the connection, make sure it satisfies Cx=b.

   
    n_constraints = length(fids_list);
    is_constrained = zeros(m.nF, 1);
    is_constrained(fids_list) = true;
    thetas = nan(m.nF, 1);
    for i = 1:n_constraints
        thetas(fids_list(i)) = thetas_list(i);
    end
    tree = cell(m.nF, 1);
    FE = m.FEAdj;
    FF = m.FFAdj;
    EF = m.EFAdj;
    
    if nargin < 4
        [~, frame_diffs] = create_local_frames(m);
    end

    FF_to_frame_diff = sparse(m.nF, m.nF);
    be = m.boundaryEdge;
    for eid = 1:m.nE
        if be(eid)
            continue;
        end
        fi = EF(eid, 1);
        fj = EF(eid, 2);
        FF_to_frame_diff(fi, fj) = frame_diffs(eid);
        FF_to_frame_diff(fj, fi) = -frame_diffs(eid);
    end

    
    %%%%%%%%%%%%%%
    % Build tree %
    %%%%%%%%%%%%%%
    
    import java.util.LinkedList
    q = LinkedList();
    f0 = fids_list(1);
    q.addLast(f0);
    visited = zeros(m.nF, 1);
    visited(f0) = true;
    tree{f0} = Nodetc(-1, -1, -1);
    
    while ~q.isEmpty()
        fi = q.removeFirst();
        
        for k = 1:3
            fj = FF(fi, k);
            eij = FE(fi, k);
            
            if fj <= 0 || visited(fj)
                continue;
            end
            
            if m.oriented_edge_is_in_face(fj, eij)
                sign = -1;
            else
                sign = 1;
            end
            
            tree{fj} = Nodetc(fi, eij, sign);
            
            visited(fj) = true;
            q.addLast(fj);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Add constrained paths %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    C = sparse(n_constraints-1, m.nE);
    b = zeros(n_constraints-1, 1);
    
    for k = 2:n_constraints
        fi = fids_list(k);
        thetai = thetas_list(k);
        % Parallel transport thetai up the tree untill you reach a
        % constrained face. Skip the root of the tree.
        while tree{fi}.p_fid > 0
            fj = tree{fi}.p_fid;
            eid = tree{fi}.p_eid;
            C(k-1, eid) = tree{fi}.sign;
            thetai = thetai + FF_to_frame_diff(fi, fj);
            
            if is_constrained(fj)
                thetaj = thetas(fj);
                break;
            end
            
            fi = fj;
        end
        b(k-1) = (thetaj - thetai);
    end
    
    

end





