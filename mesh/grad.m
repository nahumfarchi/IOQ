function [gradop, E] = grad(m)
    % GRAD Return the discrete gradient operator for the given mesh.
    % To calculate the gradient of a function f defined on the vertices
    % of m, use gradop*f.
    % E is the matrix defined here: 
    % http://webcourse.cs.technion.ac.il/236329/Winter2016-2017/ho/WCFiles/discrete_ops.pdf   
    
    % For each face, calculate the edge from the first vertex to the
    % second, second to thrid, and third to first.
    E12 = (m.V(m.F(:,2),:)-m.V(m.F(:,1),:))';
    E23 = (m.V(m.F(:,3),:)-m.V(m.F(:,2),:))';
    E31 = (m.V(m.F(:,1),:)-m.V(m.F(:,3),:))';
    
    % Build pi/2 in plane rotations matrix
    % R is a (3nf,3nf) matrix with (3,3) rotation matrices on the diagonal.
    % so, for example, R(1:3,1:3) rotates around the normal of face 1
    % R(4:6,4:6) rotates around the normal of face 2,
    % etc
    % Each rotation matrix is constructed using
    %     cos(theta)*I+sin(theta)[u]_x+u'*u
    % where theta is the angle and u is the axis of rotation and
    % [u]_x=[0 -u_z, u_y; u_z, 0, -u_x; -u_y, u_x, 0].
    % In our case, theta=pi/2, so this simplifies to
    %     [u]_x+u'*u
    I = 1:3*m.nF;
    Jr = repmat(1:3:3*m.nF, 3, 1);
    Jr = Jr(:);
    Fn = m.FNormals';
    % Arange axis vectors on the diagonal as columns vectors
    uu = sparse(I, Jr, Fn(:), 3*m.nF, 3*m.nF);
    uu = uu*uu';
    Fn = Fn';
    Ix = 3:3:3*m.nF;
    Jx = 2:3:3*m.nF;
    ux = sparse(Ix, Jx, Fn(:,1), 3*m.nF, 3*m.nF) + ...
         sparse(Jx, Ix, -Fn(:,1), 3*m.nF, 3*m.nF) + ...
         sparse(Jx-1, Ix, Fn(:,2), 3*m.nF, 3*m.nF) + ...
         sparse(Ix, Jx-1, -Fn(:,2), 3*m.nF, 3*m.nF) + ...
         sparse(Jx, Jx-1, Fn(:,3), 3*m.nF, 3*m.nF) + ...
         sparse(Jx-1, Jx, -Fn(:,3), 3*m.nF, 3*m.nF);
    R = ux + uu;
    
    Je = repmat(1:m.nF, 3, 1);
    Je = Je(:)';
    e12 = R*sparse(I, Je, E12(:), 3*m.nF, m.nF);
    e23 = R*sparse(I, Je, E23(:), 3*m.nF, m.nF);
    e31 = R*sparse(I, Je, E31(:), 3*m.nF, m.nF);
    
    ind = sub2ind(size(e12), I, Je);
    e12 = full(e12(ind));
    e12 = reshape(e12, 3, m.nF)';
    e23 = full(e23(ind));
    e23 = reshape(e23, 3, m.nF)';
    e31 = full(e31(ind));
    e31 = reshape(e31, 3, m.nF)';
    
    I = repmat(1:3*m.nF, 1, 3);
    J = [m.F(:,3); m.F(:,3); m.F(:,3); ...
         m.F(:,1); m.F(:,1); m.F(:,1); ...
         m.F(:,2); m.F(:,2); m.F(:,2)];
    E = sparse(I, J, [e12(:); e23(:); e31(:)], 3*m.nF, m.nV);
    
    gradop = 0.5*m.Gf_inv*E;
end

