function [Kg] = get_gaussian_curvature(mesh)
%function Gk = get_gaussian_curvature(mesh) 
%
% Computes the Gaussian curvature for each vertex as the discretization: 
%   2pi - sum of angles
%
% Input:
%   mesh object
%
% Output:
%   Gk - V x 1 vector

    %E12 = normalize_rows(mesh.V(mesh.F(:,2),:)-mesh.V(mesh.F(:,1),:));
    %E13 = normalize_rows(mesh.V(mesh.F(:,3),:)-mesh.V(mesh.F(:,1),:));
    %E21 = normalize_rows(mesh.V(mesh.F(:,1),:)-mesh.V(mesh.F(:,2),:));
    %E23 = normalize_rows(mesh.V(mesh.F(:,3),:)-mesh.V(mesh.F(:,2),:));
    %E31 = normalize_rows(mesh.V(mesh.F(:,1),:)-mesh.V(mesh.F(:,3),:));
    %E32 = normalize_rows(mesh.V(mesh.F(:,2),:)-mesh.V(mesh.F(:,3),:));
    
    %theta1 = acos(dot(E12, E13, 2));
    %theta2 = acos(dot(E23, E21, 2));
    %theta3 = acos(dot(E31, E32, 2)); 
    
    E12 = mesh.V(mesh.F(:,2),:)-mesh.V(mesh.F(:,1),:);
    E13 = mesh.V(mesh.F(:,3),:)-mesh.V(mesh.F(:,1),:);
    E21 = mesh.V(mesh.F(:,1),:)-mesh.V(mesh.F(:,2),:);
    E23 = mesh.V(mesh.F(:,3),:)-mesh.V(mesh.F(:,2),:);
    E31 = mesh.V(mesh.F(:,1),:)-mesh.V(mesh.F(:,3),:);
    E32 = mesh.V(mesh.F(:,2),:)-mesh.V(mesh.F(:,3),:);
    
    theta1 = atan2(row_norm(cross(E12, E13)), dot(E12, E13, 2));
    theta2 = atan2(row_norm(cross(E23, E21)), dot(E23, E21, 2));
    theta3 = atan2(row_norm(cross(E31, E32)), dot(E31, E32, 2));

    I1 = mesh.F(:,1);
    I2 = mesh.F(:,2);
    I3 = mesh.F(:,3);

    Kg = sparse([I1;I2;I3], ...
                [I2;I3;I1], ... 
                [theta1;theta2;theta3], ...
                mesh.nV, ...
                mesh.nV);

    Kg = (2*pi-full(sum(Kg,2)));
end

