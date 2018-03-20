function [Aeq, beq] = create_PU_system(m, wrapped, N)
%CREATE_PU_SYSTEM
nE_dual = m.niE;
nV_dual = m.nF;
n_rows = nE_dual;
n_cols = 2*nE_dual + nV_dual;
offset = 2*nE_dual;

EF = m.EFAdj;

I = [];
J = [];
V = [];
cost = zeros(n_cols, 1);
lb = zeros(n_cols, 1);
ub = zeros(n_cols, 1);
beq = zeros(n_rows, 1);
row = 1;
for eid = 1:m.nE
    if m.isBoundaryEdge(eid)
        continue;
    end

    fid1 = min(EF(eid, 1), EF(eid, 2));
    fid2 = max(EF(eid, 1), EF(eid, 2));

    % n_i
    I(end+1) = row;
    J(end+1) = offset + fid1;
    V(end+1) = 1;
    % n_j
    I(end+1) = row;
    J(end+1) = offset + fid2;
    V(end+1) = -1;
    % P_ij
    col = 2*(eid-1) + 1;
    I(end+1) = row;
    J(end+1) = col;
    V(end+1) = 1;
    cost(col) = 1;
    lb(col) = 0;
    ub(col) = 10;
    % Q_ij
    col = 2*(eid-1) + 2;
    I(end+1) = row;
    J(end+1) = col;
    V(end+1) = -1;
    cost(col) = 1;
    lb(col) = 0;
    ub(col) = 10;      
    % RHS
    t1 = wrapped(fid1);
    t2 = wrapped(fid2);
    kij = self.k(eid);
    period_ij = self.periods(eid);
    beq(row) = round((t2-t1+kij)/(2*pi/N), 0);
    
    row = row + 1;
end

cost(offset+1:end) = 0;
lb(offset+1:end) = -50;
ub(offset+1:end) = 50;

Aeq = full(sparse(I, J, V, n_rows, n_cols));

end

