function [F, n_flipped] = flip_edges(V, F, n_flips, angle_tol, verbose)
    % Turn
    %      2
    %   3  |  4
    %      1
    % into
    %      2
    %   3  -  4
    %      1
    % That is, (1, 2, 3) and (4, 2, 1) turn into (2, 3, 4) and (1, 4, 3).
    
    if nargin < 4
        angle_tol = pi/8;
    end
    if nargin < 5
        verbose = false;
    end
    
    n_flipped = 0;
    %for i = 1:n_flips
    while n_flipped < n_flips
        if mod(n_flipped, 100) == 0
            fprintf('n_flipped = %d\n', n_flipped);
        end
        fid1 = randi(size(F, 1));
        fid2 = find_opp_face(F, fid1);

        f1 = F(fid1, :);
        f2 = F(fid2, :);

        v1 = f1(1); 
        v2 = f1(2); 
        v3 = f1(3);
        v4 = setdiff(f2, [v1, v2]);
        
        f1_flipped = [v2, v3, v4];
        f2_flipped = [v1, v4, v3];
        
        [s1, C1] = ori(V, f1);
        s1_flipped = ori(V, f1_flipped);
        [s2, C2] = ori(V, f2);
        s2_flipped = ori(V, f2_flipped);
        
        %if abs(dot(C1, C2) - 1) > 0.1
        %    disp('skipping (non-planar faces)')
        %    continue
        %end
        if abs(acos(dot(C1,C2))) > angle_tol
            if verbose, disp('skipping (non-planar faces)'); end
            continue
        end
        
        if (s1 ~= s1_flipped || s2 ~= s2_flipped)
            if verbose, disp('skipping (inconsistent orientation)'); end
            continue
        end

        F(fid1, :) = f1_flipped;
        F(fid2, :) = f2_flipped;
        
        n_flipped = n_flipped + 1;
        
        %m = Mesh(V, F);
        %H = m.H;
    end
    
    fprintf('Flipped %d edges\n', n_flipped);
end

function fid_opp = find_opp_face(F, fid)
    f = F(fid, :);
    e = [f(1), f(2)];
    for i = 1:size(F, 1)
        if i == fid
            continue
        end
        fi = F(i, :);
        if length(intersect(fi, e)) == 2
            fid_opp = i;
            return
        end
    end
    fid_opp = -1;
end

function [y, C] = ori(V, f)
    p1 = V(f(1), :);
    p2 = V(f(2), :);
    p3 = V(f(3), :);
    v12 = p2 - p1;
    v13 = p3 - p1;
    C = cross(v12, v13);
    y = sign(C(3));
    C = C ./ norm(C);
end

% function [F] = flip_edges(F, n_flips)
%     %FLIP_EDGE Summary of this function goes here
%     %   Detailed explanation goes here
%     
%     %nv = m.nV; nf = m.nF; ne = m.nE;
%     %V = m_in.V; F = m_in.F;
%     %EF = m.EFAdj;
%     
%     %eid = randi(ne);
%     %fid1 = EF(eid, 1);
%     %fid2 = EF(eid, 2);
%     %f1 = F(fid1, :);
%     %f2 = F(fid2, :);
%     
%     
%     for i = 1:n_flips
%         fid1 = randi(size(F, 1));
%         fid2 = find_opp_face(F, fid1);
%         assert(fid1 ~= fid2)
%         if fid1 == 1642
%             disp(fid1)
%         end
%         if fid2 == 6885
%             disp(fid2)
%         end
%         if i == 131
%             disp(i)
%         end
%         f1 = F(fid1, :);
%         f2 = F(fid2, :);
%         
%         common_edge = intersect(f1, f2);
%         if ~directed_edge_is_in_face(f1, common_edge(1), common_edge(2))
%             common_edge = common_edge(end:-1:1);
%         end
%         opposite_edge = setdiff(union(f1, f2), common_edge);
%         v1_opp = setdiff(f1, common_edge);
%         v2_opp = setdiff(f2, common_edge);
%         %i1 = find(f1 == v1_opp);
%         %i2 = find(f2 == v2_opp);
%         %i2_next = mod(i2
%         %f1_flipped = [v1_opp, v2_opp, 
% 
%         f1_flipped = union(opposite_edge, common_edge(2));
%         f2_flipped = union(opposite_edge, common_edge(1));
% 
%         if directed_edge_is_in_face(f1_flipped, v2_opp, v1_opp)
%             f1_flipped = f1_flipped(end:-1:1);
%         end
%         if directed_edge_is_in_face(f2_flipped, v1_opp, v2_opp)
%             f2_flipped = f2_flipped(end:-1:1);
%         end
%         
%         %assert(orientation(V, f1) == orientation(V, f1_flipped));
%         %assert(orientation(V, f2) == orientation(V, f2_flipped));
%         assert(isempty(setdiff(union(f1, f2), union(f1_flipped, f2_flipped))))
%         assert(isempty(setdiff(union(f1_flipped, f2_flipped), union(f1, f2))))
%         assert(~isempty(setdiff(f1_flipped, f2_flipped)))
%         assert(~isempty(setdiff(f2_flipped, f1_flipped)))
%         
%         inds = setdiff(1:size(F, 1), [fid1, fid2]);
%         F = F(inds, :);
%         %assert(~face_exists(F, [3102, 3128, 3153]))
%         %assert(~face_exists(F, f1_flipped));
%         %assert(~face_exists(F, f2_flipped));
%         if ~face_exists(F, f1_flipped)
%             F = [F; f1_flipped];
%         else
%             warning('Duplicate face')
%         end
%         if ~face_exists(F, f2_flipped)
%             F = [F; f2_flipped];
%         else
%             warning('Duplicate face')
%         end
%         %F = [F; f1_flipped; f2_flipped];
%     end
%     
%     %inds = randi(ne, [n_flips, 1]);
%     %for eid = inds
%     %    fid1 = EF(e, 1);
%     %    fid2 = EF(e, 2);
%     %    f1 = F(fid1, :);
%     %    f2 = F(fid2, :);
%     %end
% end
% 
% function flip_fid = find_opp_face(F, fid)
%     v1 = F(fid, 1);
%     v2 = F(fid, 2);
%     for i = 1:size(F, 1)
%         if directed_edge_is_in_face(F(i, :), v2, v1)
%             flip_fid = i;
%             return;
%         end
%     end
%     error('?')
% end
% 
% function res = directed_edge_is_in_face(f, v1, v2)
%     if ( v1 == f(1) && v2 == f(2) ) || ...
%        ( v1 == f(2) && v2 == f(3) ) || ...
%        ( v1 == f(3) && v2 == f(1) )
%         res = true;
%     else
%         res = false;
%     end
% end
% 
% function o = orientation(V, f)
%     p1 = V(f(1), :); p2 = V(f(2), :); p3 = V(f(3), :);
%     a = p2 - p1;
%     b = p3 - p1;
%     C = cross(a, b);
%     o = sign(C(3));
% end
% 
% function y = face_exists(F, f)
%     y = false;
%     for i = 1:size(F, 1)
%         if isempty(setdiff(F(i, :), f))
%             y = true;
%             return
%         end
%     end
% end
