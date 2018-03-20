classdef MESHQ < handle
    
    properties (Access='public')
        
        name
        
        vertices
        quads
        
        nv              % #vertices
        nf              % #faces
        
        AM,valence
        
%         va              % vertex areas
%         ta              % face areas
    end
    
    methods
        
        function [ mesh ] = MESHQ( name )
                        
            mesh.name = name;
            
            [mesh.vertices, mesh.quads] = read_obj([name '.obj']);
            mesh.vertices = mesh.vertices';
            mesh.quads = mesh.quads';
                        
            mesh.nv = size(mesh.vertices,1);
            mesh.nf = size(mesh.quads,1);
            
            mesh.AM = adj_mat( mesh );
            mesh.valence = sum(mesh.AM,2);
        end
        
        function [ AM ] = adj_mat( mesh )
            T = double( mesh.quads );
            
            I = [T(:,2);T(:,3);T(:,4);T(:,1)];
            J = [T(:,3);T(:,4);T(:,1);T(:,2)];
            S = ones(size(I(:),1),1);
            AM = sparse(I,J,S,mesh.nv,mesh.nv);
        end
        
%         function [ ta ] = face_areas( mesh )
%             X = mesh.vertices;
%             T = mesh.triangles;
%             
%             P1 = X(T(:,1),:) - X(T(:,2),:);
%             P2 = X(T(:,1),:) - X(T(:,3),:);
%             
%             ta = mesh.normv( cross( P1, P2 ) ) / 2;
%         end
%         
%         function [ va ] = vertex_areas( mesh )
%             va = full( sum( mass_matrix(mesh), 2 ));
%         end
    end
    
    methods (Static)
    end
end