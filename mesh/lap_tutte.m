function [lapop] = lap_tutte(m)
    %W = (m.VVAdj' .* sparse(1:m.nV, 1:m.nV, 1./m.getValence(), m.nV, m.nV))';
    %lapop = spdiags(ones(m.nV, 1), [0], m.nV, m.nV) - W;
    lapop = spdiags(m.getValence(), [0], m.nV, m.nV) - m.VVAdj;
end

