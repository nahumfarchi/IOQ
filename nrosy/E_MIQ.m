function [E, E_edges] = E_MIQ(m, thetas, r, periods, N)
    %function [E] = E_MIQ(m, thetas, r, periods, N)
    %
    % Calculate the MIQ energy: 
    %   sum (thetas(i) + r(e_ij) + 2*pi*periods(e_ij)/N - thetas(j))^2
    %
    % Input:
    %   m - mesh object
    %   thetas - angle per face
    %   r - angle differences between frames of adjacent faces
    %   periods - period jumps between adjacent faces
    %   N - field degree
    %
    % Output:
    %   E - MIQ energy

    E_edges = zeros(m.nE, 1);
    for eid = 1:size(m.EVAdj,1)
        if m.isBoundaryEdge(eid)
            continue
        end
        i = m.EFAdj(eid, 1);
        j = m.EFAdj(eid, 2);
        E_edges(eid) = (thetas(i) + r(eid) + periods(eid)*2*pi/N - thetas(j))^2;
        %E = 
        %E = E + (thetas(i) + r(eid) + periods(eid)*2*pi/N - thetas(j))^2;
    end
    
    E = sum(E_edges);

end

