function [L] = lap_cot(mesh)
    %LAP_COT Return the cotan Laplacian matrix of the given mesh.
    %vf_1ring = mesh.vf_1ring;
    L = sparse(mesh.nV, mesh.nV);
    
    E12 = normalize_rows(mesh.V(mesh.F(:,2),:)-mesh.V(mesh.F(:,1),:));
    E13 = normalize_rows(mesh.V(mesh.F(:,3),:)-mesh.V(mesh.F(:,1),:));
    E21 = normalize_rows(mesh.V(mesh.F(:,1),:)-mesh.V(mesh.F(:,2),:));
    E23 = normalize_rows(mesh.V(mesh.F(:,3),:)-mesh.V(mesh.F(:,2),:));
    E31 = normalize_rows(mesh.V(mesh.F(:,1),:)-mesh.V(mesh.F(:,3),:));
    E32 = normalize_rows(mesh.V(mesh.F(:,2),:)-mesh.V(mesh.F(:,3),:));
    
    theta1 = acos(dot(E12, E13, 2));
    theta2 = acos(dot(E23, E21, 2));
    theta3 = acos(dot(E31, E32, 2));
    
    cot1 = cot(theta1);
    cot2 = cot(theta2);
    cot3 = cot(theta3);
    
%     E12 = mesh.V(mesh.F(:,2),:)-mesh.V(mesh.F(:,1),:);
%     E13 = mesh.V(mesh.F(:,3),:)-mesh.V(mesh.F(:,1),:);
%     E21 = mesh.V(mesh.F(:,1),:)-mesh.V(mesh.F(:,2),:);
%     E23 = mesh.V(mesh.F(:,3),:)-mesh.V(mesh.F(:,2),:);
%     E31 = mesh.V(mesh.F(:,1),:)-mesh.V(mesh.F(:,3),:);
%     E32 = mesh.V(mesh.F(:,2),:)-mesh.V(mesh.F(:,3),:);
%     
%     a_squared = sum(E12.^2, 2);
%     b_squared = sum(E13.^2, 2);
%     c_squared = sum(E23.^2, 2);
%     a = sqrt(a_squared);
%     b = sqrt(b_squared);
%     c = sqrt(c_squared);
%     
%     cos1 = 0.5 * (c_squared - a_squared - b_squared) ./ (a .* b);
%     sin1 = 2 * mesh.AF ./ (a .* b);
%     cot1 = cos1 ./ sin1;
%     
%     cos2 = 0.5 * (b_squared - c_squared - a_squared) ./ (a .* c);
%     sin2 = 2 * mesh.AF ./ (a .* c);
%     cot2 = cos2 ./ sin2;
%     
%     cos3 = 0.5 * (a_squared - b_squared - c_squared) ./ (b .* c);
%     sin3 = 2 * mesh.AF ./ (b .* c);
%     cot3 = cos3 ./ sin3;
    
    assert(~any(isnan(cot1)))
    assert(~any(isnan(cot2)))
    assert(~any(isnan(cot3)))
    
    for f = 1:mesh.nF
       src = mesh.F(f, 2);
       dst = mesh.F(f, 3);
       L(src, dst) = L(src, dst) - cot1(f);
       L(dst, src) = L(dst, src) - cot1(f);
       L(src, src) = L(src, src) + cot1(f);
       L(dst, dst) = L(dst, dst) + cot1(f);
       
       src = mesh.F(f, 1);
       dst = mesh.F(f, 3);
       L(src, dst) = L(src, dst) - cot2(f);
       L(dst, src) = L(dst, src) - cot2(f);
       L(src, src) = L(src, src) + cot2(f);
       L(dst, dst) = L(dst, dst) + cot2(f);
       
       src = mesh.F(f, 1);
       dst = mesh.F(f, 2);
       L(src, dst) = L(src, dst) - cot3(f);
       L(dst, src) = L(dst, src) - cot3(f);
       L(src, src) = L(src, src) + cot3(f);
       L(dst, dst) = L(dst, dst) + cot3(f);
    end
    
    L = 0.5 * L;

end

