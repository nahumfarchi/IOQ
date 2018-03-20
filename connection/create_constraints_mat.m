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
    F = m.F; nf = m.nF; nv = m.nV; ne = m.nE;
    
    if nargin < 4
        [~, frame_diffs] = create_local_frames(m);
    end

    FF_to_frame_diff = sparse(nf, nf);
    Eds = sparse(nf, nf);
    be = m.boundaryEdge;
    for eid = 1:m.nE
        if be(eid)
            continue;
        end
        fi = EF(eid, 1);
        fj = EF(eid, 2);
        FF_to_frame_diff(fi, fj) = frame_diffs(eid);
        FF_to_frame_diff(fj, fi) = -frame_diffs(eid);
        Eds(fi, fj) = eid;
        Eds(fj, fi) = -eid;
    end

    
    %%%%%%%%%%%%%%
    % Build tree %
    %%%%%%%%%%%%%%
    
%     import java.util.LinkedList
%     q = LinkedList();
%     f0 = fids_list(1);
%     q.addLast(f0);
%     visited = zeros(m.nF, 1);
%     visited(f0) = true;
%     tree{f0} = Nodetc(-1, -1, -1);
%     
%     while ~q.isEmpty()
%         fi = q.removeFirst();
%         
%         for k = 1:3
%             fj = FF(fi, k);
%             eij = FE(fi, k);
%             
%             if fj <= 0 || visited(fj)
%                 continue;
%             end
%             
%             if m.oriented_edge_is_in_face(fj, eij)
%                 sign = -1;
%             else
%                 sign = 1;
%             end
%             
%             tree{fj} = Nodetc(fi, eij, sign);
%             
%             visited(fj) = true;
%             q.addLast(fj);
%         end
%     end

    % Primal graph
    s = reshape(F,[],1);
    t = reshape(F(:,[2,3,1]),[],1);
    Gp = digraph(s,t); Ap = adjacency(Gp); Gp = graph(Ap);

    ne2 = sum(sum(Ap));
    [ii,jj] = find(Ap);
    Ep = sparse(ii,jj,1:ne2,nv,nv);

    % Dual graph
    F1 = sparse(s,t,[1:nf,1:nf,1:nf]',nv,nv);
    F2 = sparse(t,s,[1:nf,1:nf,1:nf]',nv,nv);
    [ii,jj,ss] = find(Ep); ll = sub2ind(size(Ep),ii,jj);
    Ed = sparse(F1(ll),F2(ll),ss,nf,nf);
    Ad = 1000*double(Ed ~= 0);
    
    for k = 1:n_constraints
        fi = fids_list(k);
        for fj = FF(fi, :)
            Ad(fi, fj) = 1;
            Ad(fj, fi) = 1;
        end
    end
   
    Gd = graph(Ad);
    
    % Find distances between the constrained faces, and use those
    % as the edge weights
%     D = distances(Gd, fids_list, fids_list);
%     assert(norm(D-D', 'fro') < 1e-10)
%     Gc = graph(D, 'upper');
%     [Tc, predc] = minspantree(Gc);   
    
    %if (sum(sum(Ad)) ~= sum(sum(Ap)))
    %    error('?');
    %end
    
    [Td, pred] = minspantree(Gd, 'method', 'sparse', 'root', fids_list(1)); 
    %Adt = adjacency(Td);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Add constrained paths %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    C = sparse(n_constraints-1, m.nE);
    b = zeros(n_constraints-1, 1);
    
%     for k = 2:n_constraints
%         fi = fids_list(k);
%         thetai = thetas_list(k);
%         % Parallel transport thetai up the tree untill you reach a
%         % constrained face. Skip the root of the tree.
%         while tree{fi}.p_fid > 0
%             fj = tree{fi}.p_fid;
%             eid = tree{fi}.p_eid;
%             C(k-1, eid) = tree{fi}.sign;
%             thetai = thetai + FF_to_frame_diff(fi, fj);
%             
%             if is_constrained(fj)
%                 thetaj = thetas(fj);
%                 break;
%             end
%             
%             fi = fj;
%         end
%         b(k-1) = (thetaj - thetai);
%     end
    
     for k = 2:n_constraints
        fi = fids_list(k);
        thetai = thetas_list(k);
        if pred(fi) == 0
            error('fi should not be root')
        end
        % parallel transport up the tree until you meet a constrained face
        while pred(fi) > 0 
            fj = pred(fi);
            eid = abs(Eds(fi, fj));
            assert(eid ~= 0, 'fi and fj are not neighbors')
            C(k-1, eid) = sign(Eds(fi, fj));
            thetai = thetai + FF_to_frame_diff(fi, fj);
            
            if is_constrained(fj)
                thetaj = thetas(fj);
                break;
            end
            
            fi = fj;
        end
        b(k-1) = -(thetaj - thetai);
    end

%     for k = 2:n_constraints
%         % Find path shortest path between two constraints on the mst
%         f_start = fids_list(k);
%         predk = predc(k);
%         f_end = fids_list(predk);
%         path = shortestpath(Gd, f_start, f_end);
%         
%         theta_end = thetas_list(predk);
%         % parallel transport theta_start to theta_end along the path
%         theta_start = thetas_list(k);
%         for j = 1:length(path)-1
%             f1 = path(j);
%             f2 = path(j+1);
%             eid = abs(Eds(f1, f2));
%             C(k-1, eid) = sign(Eds(f1, f2));
%             theta_start = theta_start + FF_to_frame_diff(f1, f2);
%         end
%         b(k-1) = theta_end - theta_start;
%     end

end





