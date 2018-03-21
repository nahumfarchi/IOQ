classdef SHAPEOP
    %SHAPEOP Summary of this class goes here
    %   Detailed explanation goes here
    % TODO: mesh,vf
        
    properties
%         Sf              % shape operator-per-face (vector)
%         sf              % shape operator-per-face (tensor)
    end
    
    methods (Static)
        
        function [ Nv ] = vertex_normals( mesh )
            
            TA = spdiags(mesh.ta,0,mesh.nf,mesh.nf);
            
            I = double( repmat(mesh.triangles(:),3,1) );
            J = repmat(1:3,3*mesh.nf,1); J = J(:);
            S = repmat(TA*mesh.Nf,3,1); S = S(:);
            
            Nv = full( sparse(I,J,S,mesh.nv,3) );
            Nv = MESH.normalize_vf( Nv );
        end
        
        function [ V, D ] = shape_operator( mesh )
            [V,D] = SHAPEOP.shape_operator_R( mesh );
        end
        
        % similar to Rusinkiewicz, 2004
        % "Estimating curvatures and their derivatives on triangle meshes"
        function [ V, D  ] = shape_operator_R( mesh )
            
            HITAR = .5 * sqrt( mesh.ITA ) * mesh.R;
            
            Nv = SHAPEOP.vertex_normals( mesh );

            gE = cell(3,1);
            gE{1} = reshape(HITAR * mesh.E1(:),[],3);
            gE{2} = reshape(HITAR * mesh.E2(:),[],3);
            gE{3} = reshape(HITAR * mesh.E3(:),[],3);

            V = zeros(mesh.nf,3);
            D = zeros(mesh.nf,1);
            for j = 1:mesh.nf
                s = 0;
                for i = 1:3
                    s = s + Nv(mesh.triangles(j,i),:)' * gE{i}(j,:);
                end
                
                % project into the tangent plane
                p = eye(3) - mesh.Nf(j,:)'*mesh.Nf(j,:);
                s = p * s;
                
                % symmetrize
                s = (s+s')/2;
                
                % compute eigenvectors
                [v,d] = eig( s ); dd = diag(d);
                
                % eliminate normal eigenvector
                e = [mesh.E1(j,:)' mesh.E2(j,:)'];
                [~,i] = min( MESH.normv( v'*e ) );
                v(:,i) = []; dd(i) = [];
                
                [~,i] = min( dd );
                V(j,:) = v(:,i)';
                D(j) = ( dd(1) - dd(2) ).^2;
            end
        end
        
        % based on Cohen-Steiner and Morvan, 2003
        function [ V, D ] = shape_operator_CSM( mesh )
            X = mesh.vertices; T = double( mesh.triangles );
            
            % adjacency matrix
            I = [T(:,2);T(:,3);T(:,1)];
            J = [T(:,3);T(:,1);T(:,2)];
            S = [1:mesh.nf,1:mesh.nf,1:mesh.nf];
            E = sparse(I,J,S,mesh.nv,mesh.nv);

            % vertex-based shape operator
            SX = zeros(9,mesh.nv);
            for i = 1:mesh.nv
                [~,n,~] = find( E(i,:) );
                sx = 0;
                for j = 1:numel(n)
                    e = X(i,:) - X(n(j),:);
                    
                    ei = mesh.v2e(min(i,n(j)),max(i,n(j)));
                    ni = mesh.Nf(abs( mesh.e2t( ei, 1 ) ),:);
                    nj = mesh.Nf(abs( mesh.e2t( ei, 2 ) ),:);
                    beta = acos( dot(ni,nj) );
                    
                    sx = sx + beta*(e'*e)/norm(e);
                end
                sx = .5*sx/mesh.va(i);
                SX(:,i) = sx(:);
            end
            
            V = zeros(mesh.nf,3);
            D = zeros(mesh.nf,1);
            
            % face-based shape operator
            for j = 1:mesh.nf
                st = reshape( sum( SX(:,T(j,:)), 2 ) / 3, [], 3 );
                
                % project into the tangent plane
                p = eye(3) - mesh.Nf(j,:)'*mesh.Nf(j,:);
                st = p * st;
                
                % symmetrize
                st = (st+st')/2;
                
                % compute eigenvectors
                [v,d] = eig( st ); v = real(v); d = real(d);
                dd = diag(d);
                
                % eliminate normal eigenvector
                e = [mesh.E1(j,:)' mesh.E2(j,:)'];
                [~,i] = min( MESH.normv( v'*e ) );
                v(:,i) = []; dd(i) = [];
                
                [~,i] = min( dd );
                V(j,:) = v(:,i)';
                D(j) = ( dd(1) - dd(2) )^2;
            end
        end
    end
    
end

