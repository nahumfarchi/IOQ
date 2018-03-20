function Ad = generator_angle_defects(mesh, verbose)
    if nargin < 2
        verbose = false;
    end
    F = mesh.F; V = mesh.V;
    cycles = mesh.generator_cycles;
    ng2 = numel(cycles);
    assert(ng2 == mesh.genus*2);
    if ng2 == 0
        Ad = [];
        return;
    end
    Ad = zeros(ng2, 1);
    for i = 1:ng2
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
end

% function [gen_Ad] = generator_angle_defects(m, local_frames, frame_diffs)
% %GENERATOR_ANGLE_DEFECTS Summary of this function goes here
% %   Detailed explanation goes here
% 
%     debug = false;
%     
%     cycles = m.generator_cycles;
%     ng = numel(cycles);
%     if ng == 0
%         gen_Ad = [];
%         return;
%     end
% 
%     if nargin < 3
%         [local_frames, frame_diffs] = create_local_frames(m);
%     end
%     
%     if debug
%         %[local_frames, frame_diffs] = create_local_frames(m);
%         scale = 0.5 * m.avg_length;
%     end
% 
%     EF = m.EFAdj;
%     FF_to_frame_diff = sparse(m.nF, m.nF);
%     be = m.boundaryEdge;
%     for eid = 1:m.nE
%         if be(eid)
%             continue;
%         end
%         fi = EF(eid, 1);
%         fj = EF(eid, 2);
%         FF_to_frame_diff(fi, fj) = frame_diffs(eid);
%         FF_to_frame_diff(fj, fi) = -frame_diffs(eid);
%     end
%     
%     gen_Ad = zeros(ng, 1);
%     for i = 1:ng
%         cy = cycles{i};
%         theta_start = 0;
%         theta_end = 0;
%         vec_start = local_angle_to_gvector(cy(1), theta_start, local_frames);
% 
%         if debug
%             figure()
%             m.draw('FaceAlpha', 0.2, 'Scale', scale)
%             hold on
%             P = sum(m.V(m.F(cy(1), :), :), 1) ./ 3;
%             arrow(P, P+0.1*vec_start, 'Length', 1, 'FaceColor', 'r')
%         end
%         
%         for j = 1:length(cy)-1
%             f1 = cy(j);
%             f2 = cy(j+1);
%             theta_end = theta_end + FF_to_frame_diff(f1, f2);
%             vec_end = local_angle_to_gvector(f2, theta_end, local_frames);
%             
%             if debug
%                 P = sum(m.V(m.F(f2, :), :), 1) ./ 3;
%                 vec = local_angle_to_gvector(f2, theta_end, local_frames);
%                 arrow(P, P+0.1*vec, 'Length', 1)
%                 if j == length(cy) - 1
%                     %vec_end = vec;
%                     a1 = vec_vec_angle(vec_end, vec_start);
%                     a2 = mod(theta_end, 2*pi) - pi;
%                     %assert(abs(theta_end - vec_vec_angle(vec_end, vec_start)) < 1e-10)
%                     assert(abs(a1 - a2) < 1e-10)
%                 end
%             end
%         end
%         gen_Ad(i) = mod(theta_end, 2*pi) - pi;
%         %gen_Ad(i) = vec_vec_angle(vec_end, vec_start);
%         
%         if debug
%             hold off
%         end
%     end
% 
% end
% 
