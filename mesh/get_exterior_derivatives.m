function [d0, d1] = get_exterior_derivatives(mesh)
    %function [d0, d1] = get_exterior_derivatives(mesh)
    %
    % Computes the exterior derivative operators on (primal) 0-forms and 1-forms of the
    % given mesh.
    %
    % Input:
    %   mesh object
    %
    % Output:
    %   d0 - E x V sparse matrix constructed as follows: for each edge (vid1, vid2) 
    %   with edge eid and vid1 < vid2, then 
    %       d0(eid, vid1) = -1 
    %       d0(eid, vid2) = +1
    %   d1 - F x E sparse matrix constructed as follows: for rach face fid with 
    %   edges eid1, eid2, eid3, then
    %       d1(fid, eid1) = +-1
    %       d2(fid, eid2) = +-1
    %       d3(fid, eid3) = +-1
    %   where the sign is decided by whether the edge agrees with the
    %   orientation of the face (+1) or not (-1).

    EV = mesh.EVAdj;
    nE = mesh.nE;
    nV = mesh.nV;

    if nargout > 0
        I = [1:nE, 1:nE];
        J = [EV(:, 1); EV(:, 2)];
        S = [-ones(nE, 1); ones(nE, 1)];
        d0 = sparse(I, J, S, nE, nV, 2*nE);
        
        %tic
        %d02 = sparse([], [], [], mesh.nE, mesh.nV, 2*mesh.nE);
        %for eid = 1:mesh.nE
        %    v1 = EV(eid, 1);
        %    v2 = EV(eid, 2);
        %    d02(eid, v1) = -1;
        %    d02(eid, v2) = 1;
        %end
        %toc
        
        %check_norm('full(d0)', 'full(d02)')
    end
    
    if nargout > 1
        % TODO remove VV_to_eids2
        VV_to_eids2 = sparse([], [], [], mesh.nV, mesh.nV, mesh.nE);
        for eid = 1:mesh.nE
            v1 = EV(eid, 1);
            v2 = EV(eid, 2);
            assert(v1 < v2);
            VV_to_eids2(v1, v2) = eid;
            VV_to_eids2(v2, v1) = eid;
        end
        
        v1 = EV(:, 1);
        v2 = EV(:, 2);
        VV_to_eids = sparse([v1; v2], [v2; v1], [1:mesh.nE, 1:mesh.nE]', mesh.nV, ...
            mesh.nV, 2*mesh.nE);
        assert(norm(VV_to_eids2 - VV_to_eids, 'fro') == 0)
        
        d1 = sparse([], [], [], mesh.nF, mesh.nE, 3*mesh.nF);
        for fid = 1:mesh.nF
            for i = 1:3
                v1 = mesh.F(fid, i);
                v2 = mesh.F(fid, mod(i,3)+1);

                sgn = 1;
                if v1 > v2
                    tmp = v1;
                    v1 = v2;
                    v2 = tmp;
                    sgn = -1;
                end

                eid = VV_to_eids(v1, v2);

                d1(fid, eid) = sgn;
            end
        end
    end
end

