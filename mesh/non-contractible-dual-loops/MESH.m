classdef MESH < handle
    
    properties (Access='public')
        
        name
        
        vertices
        triangles
        edges
        
        e2t, t2e, v2e
        
        nv              % #vertices
        nf              % #faces
        ne              % #edges
        
        va              % vertex areas
        ta              % face areas
        ea              % edge areas
        
        VA, IVA         % (lumped) vertices mass matrix and inverse
        TA, ITA         % (lumped) triangles mass matrix and inverse
        EA              % (lumped) edges mass matrix (hmm, inverse?)
        
        Nf              % normal-per-face
        Nv              % normal-per-vertex
        R               % in-plane rotation operator wrt normal
                
        G               % gradient of vertex functions
        D               % divergence of face fields
        C               % curl of face fields
        GE              % gradient of edge functions
        DE              % divergence of face fields
        
        Iv2f            % interpolate vertices to faces
        If2v            % interpolate faces to vertices        
        Ie2f            % interpolate edges to faces
        If2e            % interpolate faces to edges
        
        E1, E2, E3      % triangle edges
        F1, F2          % frame per face
        EB, EBI         % an edge based basis for vfs.
        
        mel             % mean edge length
    end
    
    properties (Access='public')
%         Sf              % shape operator-per-face (vector)
%         Sf2             % second shape operator-per-face (vector)
%         sf              % shape operator-per-face (tensor)
%         sf2             % second shape operator-per-face (tensor)
%         Hv              % vertex mean curvature
%         Hf              % face mean curvature
%         Kv              % vertex Gaussian curvature
%         Kf              % face Gaussian curvature
%         
%         Tv              % vertex total curvature
%         HS              % mean curvature + shape operator (tensor)
%         z               % up direction
%         cos_theta       % Nv \cdot z
    end
    
    methods
        function [ mesh ] = MESH( meshname, v, t)
            if nargin < 1; meshname = 'sphere_s3'; end
            mesh.name = meshname;

            if nargin < 3   
                [v, t] = MESH_IO.roff([meshname '.off']);
            end
            
            mesh.vertices = v;
            mesh.triangles = t;
            
            % mesh data structures
                        
            mesh.nv = size(mesh.vertices,1);
            mesh.nf = size(mesh.triangles,1);
            
            [mesh.edges,mesh.e2t,mesh.t2e,mesh.v2e] = nc_data( mesh );
            mesh.ne = size(mesh.edges,1);
            
            % simplices' areas and mass matrices
            mesh.ta = face_areas( mesh );
            mesh.va = vertex_areas( mesh );
            mesh.ea = edge_areas( mesh );
            [mesh.VA,mesh.IVA,mesh.TA,mesh.ITA,mesh.EA] = mass_matrices( mesh );
            
            % simplices' normals
            mesh.Nf = face_normals( mesh );
            mesh.Nv = vertex_normals( mesh );
            
            % geometric and differential operators
            [mesh.E1,mesh.E2,mesh.E3] = face_edges( mesh );
            mesh.R = rot( mesh );
                        
            mesh.G = grad( mesh );
            mesh.D = div( mesh );
            
            mesh.GE = grade( mesh );
            mesh.DE = dive( mesh );
            
            mesh.C = curl( mesh );
            
            % interpolants
            mesh.If2v = f2v( mesh );
            mesh.Iv2f = v2f( mesh );
            
            mesh.If2e = f2e( mesh );
            mesh.Ie2f = e2f( mesh );
            
            [mesh.F1, mesh.F2, mesh.EB, mesh.EBI] = eb( mesh );
        end
        
        function [ edges, e2t, t2e, v2e ] = nc_data( mesh )
            T = double( mesh.triangles );
            
            I = [T(:,2);T(:,3);T(:,1)];
            J = [T(:,3);T(:,1);T(:,2)];
            S = [1:mesh.nf,1:mesh.nf,1:mesh.nf];
            E = sparse(I,J,S,mesh.nv,mesh.nv);
            
            Elisto = [I,J];
            sElist = sort(Elisto,2);
            s = (MESH.normv(Elisto - sElist) > 1e-12);
            t = S'.*(-1).^s;
            [edges,une] = unique(sElist, 'rows');
            e2t = zeros(length(edges),4);
            t2e = zeros(mesh.nf,3);
            for m=1:length(edges)
                i = edges(m,1); j = edges(m,2);
                t1 = t(une(m));
                t2 = -(E(i,j) + E(j,i) - abs(t1))*sign(t1);
                e2t(m,1:2) = [t1, t2];
                f = T(abs(t1),:); loc = find(f == (sum(f) - i - j));
                t2e(abs(t1),loc) = m*sign(t1);
                e2t(m,3) = loc;
                if t2 ~= 0
                    f = T(abs(t2),:); loc = find(f == (sum(f) - i - j));
                    t2e(abs(t2),loc) = m*sign(t2);
                    e2t(m,4) = loc;
                end
            end

            v2e = sparse(edges(:,1),edges(:,2),1:length(edges),mesh.nv,mesh.nv);
        end
        
        function [ ta ] = face_areas( mesh )
            X = mesh.vertices;
            T = mesh.triangles;
            
            P1 = X(T(:,1),:) - X(T(:,2),:);
            P2 = X(T(:,1),:) - X(T(:,3),:);
            
            ta = mesh.normv( cross( P1, P2 ) ) / 2;
        end
        
        function [ va ] = vertex_areas( mesh )
            va = full( sum( mass_matrix(mesh), 2 ));
        end
        
        function [ M ] = mass_matrix( mesh )
            T = double( mesh.triangles ); 
            
            I = [T(:,1);T(:,2);T(:,3)];
            J = [T(:,2);T(:,3);T(:,1)];
            Mij = 1/12*[mesh.ta; mesh.ta; mesh.ta];
            Mji = Mij;
            Mii = 1/6*[mesh.ta; mesh.ta; mesh.ta];
            In = [I;J;I];
            Jn = [J;I;I];
            Mn = [Mij;Mji;Mii];
            M = sparse(In,Jn,Mn,mesh.nv,mesh.nv);
            
% X = mesh.vertices; T = mesh.triangles;
% nv = mesh.nv;
% 
% L1 = normv(X(T(:,2),:)-X(T(:,3),:));
% L2 = normv(X(T(:,1),:)-X(T(:,3),:));
% L3 = normv(X(T(:,1),:)-X(T(:,2),:));
% A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2.*L2.*L3);
% A2 = (L1.^2 + L3.^2 - L2.^2) ./ (2.*L1.*L3);
% A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2.*L1.*L2);
% A = [A1,A2,A3];
% A = acos(A);
% EL = [L1, L2, L3];
% 
% ar = (1/8)*cot(A).*EL.^2;
% I = [T(:,1);T(:,2);T(:,3)];
% J = [T(:,2);T(:,3);T(:,1)];
% S = [ar(:,3);ar(:,1);ar(:,2)];
% In = [I;J];
% Jn = In;
% Sn = [S;S];
% M = sparse(In,Jn,Sn,nv,nv); 
            
        end
        
        function [ ea ] = edge_areas( mesh )
            T = double( mesh.triangles );

            I = [T(:,1);T(:,2);T(:,3)];
            J = [T(:,2);T(:,3);T(:,1)];
            S = [1:mesh.nf,1:mesh.nf,1:mesh.nf]';
            E = sparse(I,J,S',mesh.nv,mesh.nv);

            RN = E( [ sub2ind( size(E), T(:,3), T(:,2) ), ...
                      sub2ind( size(E), T(:,1), T(:,3) ), ...
                      sub2ind( size(E), T(:,2), T(:,1) ) ] );
            RN = full( RN );

            I = repmat((1:mesh.nf),3,1); I = I(:);
            J = RN'; J = J(:);
            S = repmat(mesh.ta/3,1,3)'; S = S(:);
            S( J == 0 ) = []; I( J == 0 ) = []; J( J == 0 ) = []; % bdry
            W = sparse(I,J,S,mesh.nf,mesh.nf);
            
            E2T = abs( mesh.e2t(:,1:2) );
            
            % boundary handling
            II = find(E2T(:,1) ~= 0 & E2T(:,2) ~= 0);
            IB = find(E2T(:,1) == 0 | E2T(:,2) == 0);
            
            i1 = sub2ind(size(W),E2T(II,1),E2T(II,2));
            i2 = sub2ind(size(W),E2T(II,2),E2T(II,1));
            ea = zeros(mesh.ne,1);
            ea(II) = full( W(i1) + W(i2) );
            ea(IB) = mesh.ta(sum(E2T(IB,:),2))/3;
        end
        
        function [ VA, IVA, TA, ITA, EA ] = mass_matrices( mesh )
            sv = mesh.nv; sf = mesh.nf; se = mesh.ne; 
            VA = spdiags(mesh.va,0,sv,sv);
            IVA = spdiags(1./mesh.va,0,sv,sv);
            TA = spdiags([mesh.ta; mesh.ta; mesh.ta],0,3*sf,3*sf);
            ITA = spdiags(1./[mesh.ta; mesh.ta; mesh.ta],0,3*sf,3*sf);
            EA = spdiags(mesh.ea,0,se,se);
        end
        
        function [ Nf ] = face_normals( mesh )
            
            X = mesh.vertices;
            T = mesh.triangles;
            
            P1 = X(T(:,1),:) - X(T(:,2),:);
            P2 = X(T(:,1),:) - X(T(:,3),:);
            
            Nf = cross( P1, P2 );
            Nf = MESH.normalize_vf( Nf );
        end
        
        function [ Nv ] = vertex_normals( mesh )
            if exist([mesh.name '_normals.off']) == 2
                [Nv, ~] = MESH_IO.roff([mesh.name '_normals.off']);
                fprintf('Read mesh normals from file\n');
                return;
            end
            
%             TA = spdiags(mesh.ta,0,mesh.nf,mesh.nf);
            sf = mesh.nf;
            
            I = double( repmat(mesh.triangles(:),3,1) );
            J = repmat(1:3,3*mesh.nf,1); J = J(:);
            
            % S = repmat(mesh.Nf,3,1); S = S(:);
            S = repmat(mesh.TA(1:sf,1:sf)*mesh.Nf,3,1); S = S(:);
            % S = repmat(mesh.ia.*mesh.Nf,3,1); S = S(:);
            
            Nv = full( sparse(I,J,S,mesh.nv,3) );
            Nv = MESH.normalize_vf( Nv );
        end
        
        function [ E1, E2, E3 ] = face_edges( mesh )
            X = mesh.vertices;
            T = mesh.triangles;
            
            E1 = X(T(:,3),:) - X(T(:,2),:);
            E2 = X(T(:,1),:) - X(T(:,3),:);
            E3 = X(T(:,2),:) - X(T(:,1),:);
            E = [E1; E2; E3];
            
            mesh.mel = mean( MESH.normv( E ) );
        end
        
        function [ R ] = rot( mesh )
            sf = mesh.nf; n = mesh.Nf;
            
            II = repmat((1:sf)',1,2); II = [II II+sf II+2*sf]'; II = II(:);
            JJ1 = (1:sf)'; JJ2 = JJ1+sf; JJ3 = JJ2+sf;
            JJ = [JJ2 JJ3 JJ1 JJ3 JJ1 JJ2]'; JJ = JJ(:);
            SS = [-n(:,3) n(:,2) n(:,3) -n(:,1) -n(:,2) n(:,1)]'; SS = SS(:);
            
            R = sparse(II,JJ,SS,3*sf,3*sf);
        end
        
        function [ G ] = grad( mesh )
            % G corresponds to eq. (3.9) in Polygon mesh processing book
            I = repmat(1:mesh.nf,3,1);
            II = [I(:); I(:)+mesh.nf; I(:)+2*mesh.nf];

            J = double( mesh.triangles' );
            JJ = [J(:); J(:); J(:)];

            RE1 = mesh.rotate_vf( mesh.E1 );
            RE2 = mesh.rotate_vf( mesh.E2 );
            RE3 = mesh.rotate_vf( mesh.E3 );

            S = [RE1(:) RE2(:) RE3(:)]'; SS = S(:);
            
            G = sparse(II,JJ,SS,3*mesh.nf,mesh.nv);
            G = .5 * mesh.ITA * G;
        end
        
        function [ D ] = div( mesh )
            % D corresponds to eq. (3.12) in Polygon mesh processing book
            D = - mesh.IVA * mesh.G' * mesh.TA;
        end
        
        function [ C ] = curl( mesh )
            C = - mesh.DE * mesh.R;
%             C = - mesh.D * mesh.R;
        end
        
        function [ GE ] = grade( mesh )
            X = mesh.vertices;
            E = mesh.edges;
            
            I = repmat(1:mesh.nf,3,1);
            II = [I(:); I(:)+mesh.nf; I(:)+2*mesh.nf];

            J = double( abs(mesh.t2e)' );
            JJ = [J(:); J(:); J(:)];
            
            f1 = repmat( sign(mesh.t2e(:,1)), 1, 3 );
            f2 = repmat( sign(mesh.t2e(:,2)), 1, 3 );
            f3 = repmat( sign(mesh.t2e(:,3)), 1, 3 );
            
            at2e = abs(mesh.t2e);
            EE1 = f1 .* ( X(E(at2e(:,1),1),:) - X(E(at2e(:,1),2),:) );
            EE2 = f2 .* ( X(E(at2e(:,2),1),:) - X(E(at2e(:,2),2),:) );
            EE3 = f3 .* ( X(E(at2e(:,3),1),:) - X(E(at2e(:,3),2),:) );
            
            RE1 = - mesh.rotate_vf( EE1 );
            RE2 = - mesh.rotate_vf( EE2 );
            RE3 = - mesh.rotate_vf( EE3 );
            
            S = [RE1(:) RE2(:) RE3(:)]'; SS = S(:);
            
            GE = sparse(II,JJ,SS,3*mesh.nf,mesh.ne);
            GE = mesh.ITA * GE;
        end
        
        function [ DE ] = dive( mesh )
            IEA = spdiags(1./mesh.ea,0,mesh.ne,mesh.ne);
            TAC = spdiags(repmat(mesh.ta,3,1),0,3*mesh.nf,3*mesh.nf);

            DE = - IEA * mesh.GE' * TAC;
        end
        
        function [ DF ] = divf( mesh )
            T = double( mesh.triangles );
            
            II = [T(:,2);T(:,3);T(:,1)];
            JJ = [T(:,3);T(:,1);T(:,2)];
            SS = [1:mesh.nf,1:mesh.nf,1:mesh.nf]';
            EE = sparse(II,JJ,SS',mesh.nv,mesh.nv);
            
%             L1 = MESH.normv( mesh.E1 );
%             L2 = MESH.normv( mesh.E2 );
%             L3 = MESH.normv( mesh.E3 );
%             
%             S1 = (L2.^2 + L3.^2 - L1.^2) ./ (2.*L2.*L3);
%             S2 = (L3.^2 + L1.^2 - L2.^2) ./ (2.*L3.*L1);
%             S3 = (L1.^2 + L2.^2 - L3.^2) ./ (2.*L1.*L2);
%             SS2 = [S1; S2; S3];
%             SS2 = 0.5*cot( acos(SS2) );
% %             SS2 = repmat(mesh.ta/3,3,1);
%             EE2 = sparse(II,JJ,SS2,mesh.nv,mesh.nv);
%             ea = EE2( [ sub2ind( size(EE), T(:,3), T(:,2) ), ...
%                        sub2ind( size(EE), T(:,1), T(:,3) ), ...
%                        sub2ind( size(EE), T(:,2), T(:,1) ) ] ) + ...
%                  EE2( [ sub2ind( size(EE), T(:,2), T(:,3) ), ...
%                        sub2ind( size(EE), T(:,3), T(:,1) ), ...
%                        sub2ind( size(EE), T(:,1), T(:,2) ) ] );
%             w = repmat(ea(:),6,1);
            
            % 1-ring of each face, ordered as the edges
            TN = EE( [ sub2ind( size(EE), T(:,3), T(:,2) ), ...
                       sub2ind( size(EE), T(:,1), T(:,3) ), ...
                       sub2ind( size(EE), T(:,2), T(:,1) ) ] );
            TN = full( TN );
            
            RE1 = mesh.rotate_vf( mesh.E1 );
            RE2 = mesh.rotate_vf( mesh.E2 );
            RE3 = mesh.rotate_vf( mesh.E3 );
            reb = [RE1; RE2; RE3];
                        
            TN = TN(:); reb = reb(:);

            w = repmat(mesh.ta/3,3,1) + mesh.ta(TN)/3;
            w = repmat(1./w,3,1); w = [w; w];
             
            I = [TN; TN; TN; repmat((1:mesh.nf)',9,1)];
            J = repmat((1:mesh.nf)',3,1); J = [J; J+mesh.nf; J+2*mesh.nf];
            J = [J; J];
            S = 1/3 * w .* [reb; reb];
%             S = 1/3 * [reb; reb];
                       
            DF = sparse(I,J,S,mesh.nf,3*mesh.nf);
        end
        
        function [ op ] = fvf( mesh, vf )
            % op = [v]_{\cdot}^T * grad \in nf x nv
            
            I = repmat(1:mesh.nf,1,3); I = I(:); 
            J = (1:mesh.nf)'; J = [J; J+mesh.nf; J+2*mesh.nf];
            V = sparse(I,J,vf,mesh.nf,3*mesh.nf);

            op = V * mesh.G;
        end

        function [ op ] = ifvf( mesh, vf )
            % op = I^F_V * [v]_{\cdot}^T * grad \in nv x nv
            
            op = mesh.If2v * mesh.fvf( vf );
        end
        
        function [ op ] = ifvfe( mesh, vf )
            % op = I^F_V * [v]_{\cdot}^T * grad \in ne x nv
            
            op = mesh.If2e * mesh.fvf( vf );
        end
        
        function [ op ] = advf( mesh, vf )
            % op = div * [v]_{\cdot} \in nv x nf
            
            I = repmat(1:mesh.nf,1,3); I = I(:); 
            J = (1:mesh.nf)'; J = [J; J+mesh.nf; J+2*mesh.nf];
            V = sparse(I,J,vf,mesh.nf,3*mesh.nf);
            
            op = mesh.D * V';
        end

        function [ op ] = advfe( mesh, vf )
            % op = div * [v]_{\cdot} \in nv x ne
            
            op = mesh.advf( vf ) * mesh.If2e'; 
        end
        
        function [ op ] = adfvf( mesh, vf )
            % op = div * [v]_{\cdot} - [div(v)] * I^F_V \in nv x nf
            
            op = mesh.advf( vf ) - ...
                 spdiags( mesh.D * vf(:), 0, mesh.nv, mesh.nv ) * mesh.If2v;
        end

        function [ op ] = adfvfe( mesh, vf )
            % op = div * [v]_{\cdot} - [div(v)] * I^F_V \in nv x ne
            
            op = mesh.advfe( vf ) - ...
                 spdiags( mesh.D * vf(:), 0, mesh.nv, mesh.nv ) * mesh.If2v * mesh.If2e';
        end
        
        function [ op ] = iadfvf( mesh, vf )
            % op = div * [v]_{\cdot} * I^V_F - [div(v)] \in nv x nv
            
            op = mesh.advf( vf ) * mesh.Iv2f - ...
                 spdiags( mesh.D * vf(:), 0, mesh.nv, mesh.nv );
        end
        
        function [ op ] = fvfbar( mesh, f )
            % op = [grad f]_{\cdot}^T \in nf x 3nf
            
            gf = mesh.G * f;
            
            I = repmat(1:mesh.nf,1,3); I = I(:); 
            J = (1:mesh.nf)'; J = [J; J+mesh.nf; J+2*mesh.nf];
            GF = sparse(I,J,gf,mesh.nf,3*mesh.nf);
            
            op = GF;
        end
        
        function [ op ] = liebr( mesh, u, v )
            % op \in nv x nv
            
            fvfu = mesh.fvf( u );
            fvfv = mesh.fvf( v );
            
            adfvfu = mesh.adfvf( u );
            adfvfv = mesh.adfvf( v );
                        
            op = adfvfu*fvfv - adfvfv*fvfu;
        end

        function [ op ] = liebre( mesh, u, v )
            % op \in nv x nv
            
            fvfu = mesh.ifvfe( u );
            fvfv = mesh.ifvfe( v );
            
            adfvfu = mesh.adfvfe( u );
            adfvfv = mesh.adfvfe( v );
                        
            op = adfvfu*fvfv - adfvfv*fvfu;
        end
        
        function [ op ] = liedrv( mesh, u )
            % op in 3nf x 3nf
            sv = mesh.nv; sf = mesh.nf;
            
            I = repmat(1:sf,1,3); I = I(:); 
            J = (1:sf)'; J = [J; J+sf; J+2*sf];
            U = sparse(I,J,u,sf,3*sf);
            
            RUB1 = spdiags(repmat(u(:,1),3,1),0,3*sf,3*sf);
            RUB2 = spdiags(repmat(u(:,2),3,1),0,3*sf,3*sf);
            RUB3 = spdiags(repmat(u(:,3),3,1),0,3*sf,3*sf);
            
            spz = spalloc(3*sf,sf,1);
            
            op1 = [ mesh.D*( [U' spz spz] - RUB1 ); ...
                    mesh.D*( [spz U' spz] - RUB2 ); ...
                    mesh.D*( [spz spz U'] - RUB3 ) ]; 
            
            DU = spdiags(-mesh.D*u(:),0,sv,sv);
            DU = DU * mesh.If2v;
            
            uv = mesh.If2v * u;

            UV1 = spdiags(uv(:,1),0,sv,sv);
            UV2 = spdiags(uv(:,2),0,sv,sv);
            UV3 = spdiags(uv(:,3),0,sv,sv);
            
            spz = spalloc(sv,sf,1);
            
            op2 = [DU spz spz; spz DU spz; spz spz DU] + ...
                  [UV1 * mesh.D; UV2 * mesh.D; UV3 * mesh.D];
            
            op = op1 + op2;
        end

        function [ op ] = liedrvi( mesh, u )
            % op in 3nf x 3nf
            sv = mesh.nv; sf = mesh.nf;
            
            Ru = mesh.R*u(:);    
            I = repmat(1:mesh.nf,1,3); I = I(:); 
            J = (1:mesh.nf)'; J = [J; J+mesh.nf; J+2*mesh.nf];
            RU = sparse(I,J,Ru,mesh.nf,3*mesh.nf);
            op1 = mesh.R*mesh.G*mesh.If2v*RU;
                
            DU = spdiags(-mesh.D*u(:),0,sv,sv);
            DU = DU * mesh.If2v;
            
            uv = mesh.If2v * u;

            UV1 = spdiags(uv(:,1),0,sv,sv);
            UV2 = spdiags(uv(:,2),0,sv,sv);
            UV3 = spdiags(uv(:,3),0,sv,sv);
            
            spz = spalloc(sv,sf,1);
            
            op2 = [DU spz spz; spz DU spz; spz spz DU] + ...
                  [UV1 * mesh.D; UV2 * mesh.D; UV3 * mesh.D];
            
            op = op1 + blkdiag(mesh.Iv2f,mesh.Iv2f,mesh.Iv2f)*op2;
        end

        function [ op ] = liedrvie( mesh, u )
            % op in 3nf x 3nf
            sv = mesh.nv; sf = mesh.nf;
            
            Ru = mesh.R*u(:);    
            I = repmat(1:mesh.nf,1,3); I = I(:); 
            J = (1:mesh.nf)'; J = [J; J+mesh.nf; J+2*mesh.nf];
            RU = sparse(I,J,Ru,mesh.nf,3*mesh.nf);
            op1 = mesh.R*mesh.GE*mesh.If2e*RU;
                
            DU = spdiags(-mesh.D*u(:),0,sv,sv);
            DU = DU * mesh.If2v;
            
            uv = mesh.If2v * u;

            UV1 = spdiags(uv(:,1),0,sv,sv);
            UV2 = spdiags(uv(:,2),0,sv,sv);
            UV3 = spdiags(uv(:,3),0,sv,sv);
            
            spz = spalloc(sv,sf,1);
            
            op2 = [DU spz spz; spz DU spz; spz spz DU] + ...
                  [UV1 * mesh.D; UV2 * mesh.D; UV3 * mesh.D];
            
            op = op1 + blkdiag(mesh.Iv2f,mesh.Iv2f,mesh.Iv2f)*op2;
        end
        
        function [ op ] = liedre( mesh, u )
            % op in 3nf x 3nf
            se = mesh.ne; sf = mesh.nf;
            
            I = repmat(1:sf,1,3); I = I(:); 
            J = (1:sf)'; J = [J; J+sf; J+2*sf];
            U = sparse(I,J,u,sf,3*sf);
            
            RUB1 = spdiags(repmat(u(:,1),3,1),0,3*sf,3*sf);
            RUB2 = spdiags(repmat(u(:,2),3,1),0,3*sf,3*sf);
            RUB3 = spdiags(repmat(u(:,3),3,1),0,3*sf,3*sf);
            
            spz = spalloc(3*sf,sf,1);
            
            op1 = [ mesh.DE*( [U' spz spz] - RUB1 ); ...
                    mesh.DE*( [spz U' spz] - RUB2 ); ...
                    mesh.DE*( [spz spz U'] - RUB3 ) ]; 
            
            DU = spdiags(-mesh.DE*u(:),0,se,se);
            DU = DU * mesh.If2e;
            
            ue = mesh.If2e * u;

            UE1 = spdiags(ue(:,1),0,se,se);
            UE2 = spdiags(ue(:,2),0,se,se);
            UE3 = spdiags(ue(:,3),0,se,se);
            
            spz = spalloc(se,sf,1);
            
            op2 = [DU spz spz; spz DU spz; spz spz DU] + ...
                  [UE1 * mesh.DE; UE2 * mesh.DE; UE3 * mesh.DE];
            
            op = op1 + op2;
        end
        
        function [ op ] = hodge( mesh )
            op = mesh.D' * mesh.VA * mesh.D + mesh.C' * mesh.EA * mesh.C;
%             op = mesh.D' * mesh.VA * mesh.D + mesh.C' * mesh.VA * mesh.C;
            op = mesh.ITA * op;
            op = mesh.EB * op * mesh.EBI;
        end
        
        function [ op ] = vdm( mesh )
            sv = mesh.nv; sf = mesh.nf;
            
            T = double( mesh.triangles );
            
            I = [T(:,1);T(:,2);T(:,3)];
            J = [T(:,2);T(:,3);T(:,1)];
            S = [1:sf,1:sf,1:sf]';
            E = sparse(I,J,S',sv,sv);
            
            RN = E( [ sub2ind( size(E), T(:,3), T(:,2) ), ...
                      sub2ind( size(E), T(:,1), T(:,3) ), ...
                      sub2ind( size(E), T(:,2), T(:,1) ) ] );
            RN = full( RN );
            
            % weight matrices
            X = mesh.vertices; T = mesh.triangles;
            X = (X(T(:,1),:)+X(T(:,2),:)+X(T(:,3),:))/3;
            
            I = repmat((1:sf),3,1); I = I(:);
            J = RN'; J = J(:);
            J( J == 0 ) = I( J == 0 ); % Boundary handling
            w = exp( - MESH.snormv(X(I,:) - X(J,:)) ); % \epsilon = 1
            W = sparse(I,J,w,sf,sf);
            
            deg = sum( W, 2 );
            DI = spdiags(1./deg,0,sf,sf);
            W1 = DI * W * DI;
            deg1 = sum( W1, 2 );
            rdeg1 = repmat(deg1,2,1);
            RD1I = spdiags(1./rdeg1,0,2*sf,2*sf);            
            
            rdeg = repmat(deg,2,1);
%             D = spdiags(rdeg,0,2*sf,2*sf);
            RDI = spdiags(1./rdeg,0,2*sf,2*sf);

            NE1 = MESH.normalize_vf( mesh.E1 );
            NE2 = rotate_vf(mesh,NE1);
            
            II = (1:sf)'; II = [II II+sf];
            II = reshape(repmat(II,1,6)',[],1);
            
            JJ = reshape(RN',[],1); JJ = [JJ JJ JJ+sf JJ+sf];
            JJ = reshape(JJ',[],1);
            
            % Boundary handling
            k = find( JJ == 0 ); k = reshape(k,2,[])';
            k = [k k+repmat([2 2],size(k,1),1)]; k = reshape(k',[],1);
            II(k) = []; JJ(k) = [];
            
            SS = [];
            for i = 1:sf
                Oit = [NE1(i,:); NE2(i,:)];
           
                for j = 1:length(RN(i,:))
                    rnij = RN(i,j);
                    
                    % Boundary edges
                    if rnij==0; continue; end 
                                        
                    Oj = [NE1(rnij,:)' NE2(rnij,:)'];

                    Oij = Oit*Oj;
                    [u,s,v] = svd(Oij);
                    u = repmat(sign(diag(s))',size(u,1),1).*u;
                    
%                     wij = 1/3;
                    wij = W(i,rnij);
                    SS = [SS; reshape(wij*u*v',[],1);];
                end
            end
            S = sparse(II,JJ,SS,2*sf,2*sf);
            S1 = RDI * S * RDI;
            
            DS =  RD1I * S1;
%             op = .5*(DS+DS');
            
            I = speye(2*sf,2*sf);
            op = I-DS;
            op = .5*(op+op');
        end
        
        function [ op ] = vdm2( mesh )
            sv = mesh.nv; sf = mesh.nf;
            
            T = double( mesh.triangles );
            
            I = [T(:,1);T(:,2);T(:,3)];
            J = [T(:,2);T(:,3);T(:,1)];
            S = [1:sf,1:sf,1:sf]';
            E = sparse(I,J,S',sv,sv);
            
            RN = E( [ sub2ind( size(E), T(:,3), T(:,2) ), ...
                      sub2ind( size(E), T(:,1), T(:,3) ), ...
                      sub2ind( size(E), T(:,2), T(:,1) ) ] );
            RN = full( RN );
            
            NE1 = MESH.normalize_vf( mesh.E1 );
            NE2 = rotate_vf(mesh,NE1);
            
            II = (1:sf)'; II = [II II+sf];
            II = reshape(repmat(II,1,6)',[],1);
            
            JJ = reshape(RN',[],1); JJ = [JJ JJ JJ+sf JJ+sf];
            JJ = reshape(JJ',[],1);
            
            % Boundary handling
            k = find( JJ == 0 ); k = reshape(k,2,[])';
            k = [k k+repmat([2 2],size(k,1),1)]; k = reshape(k',[],1);
            II(k) = []; JJ(k) = [];
            
            SS = [];
            for i = 1:sf
                Oit = [NE1(i,:); NE2(i,:)];
           
                for j = 1:length(RN(i,:))
                    rnij = RN(i,j);
                    
                    % Boundary edges
                    if rnij==0; continue; end 
                                        
                    Oj = [NE1(rnij,:)' NE2(rnij,:)'];

                    Oij = Oit*Oj;
                    [u,s,v] = svd(Oij);
                    u = repmat(sign(diag(s))',size(u,1),1).*u;
                    
                    wij = 1/3;
                    SS = [SS; reshape(wij*u*v',[],1);];
                end
            end
            DS = sparse(II,JJ,SS,2*sf,2*sf);
%             op = .5*(DS+DS');
            
            I = speye(2*sf,2*sf);
            op = I-DS;
            op = .5*(op+op');
        end
        function [K] = gaussian_curv(mesh)
            X = mesh.vertices; T = mesh.triangles;
            
            % Find orig edge lengths and angles
            L1 = normv(X(T(:,2),:)-X(T(:,3),:));
            L2 = normv(X(T(:,1),:)-X(T(:,3),:));
            L3 = normv(X(T(:,1),:)-X(T(:,2),:));
            A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2.*L2.*L3);
            A2 = (L1.^2 + L3.^2 - L2.^2) ./ (2.*L1.*L3);
            A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2.*L1.*L2);
            A = [A1,A2,A3];
            A = acos(A);
            
            % Find original curvature - Korig
            I = reshape(double(T),mesh.nf*3,1);
            J = ones(mesh.nf*3,1);
            S = reshape(A,mesh.nf*3,1);
            SA = sparse(I,J,S,mesh.nv,1);    
            SA = full(SA);
            K = 2*pi*ones(mesh.nv,1) - SA;
        end
        
        function [ op ] = v2f( mesh )
            
%             ITA = spdiags(1./mesh.ta,0,mesh.nf,mesh.nf);
            sf = mesh.nf;
            op = mesh.ITA(1:sf,1:sf) * mesh.If2v' * mesh.VA;
        end
    
        function [ op ] = f2v( mesh )

            I = double( mesh.triangles ); I = I(:);
            J = repmat(1:mesh.nf,3,1)'; J = J(:);
            S = repmat(mesh.ta/3,3,1);
                        
            op = mesh.IVA * sparse(I,J,S,mesh.nv,mesh.nf);
        end
        
        function [ op ] = e2f( mesh )
%             ITA = spdiags(1./mesh.ta,0,mesh.nf,mesh.nf);
            sf = mesh.nf;
            op = mesh.ITA(1:sf,1:sf) * mesh.If2e' * mesh.EA;
        end
        
        function [ op ] = f2e( mesh )
            IEA = spdiags(1./mesh.ea,0,mesh.ne,mesh.ne);
            
            I = double( abs(mesh.t2e) ); I = I(:);
            J = repmat(1:mesh.nf,3,1)'; J = J(:);
            S = repmat(mesh.ta/3,3,1);
                        
            op = IEA*sparse(I,J,S,mesh.ne,mesh.nf);
        end
        
        function [ rvf ] = rotate_vf( mesh, vf )
            vf = reshape(vf,mesh.nf,3);
            rvf = cross( mesh.Nf, vf );
        end
        
        function [ pvf ] = project_vf( mesh, vf )
            pvf = vf - mesh.Nf .* repmat(dot(mesh.Nf,vf,2),1,3);
        end
        
        function [ F1, F2, EB, EBI ] = eb( mesh )
            
            F1 = MESH.normalize_vf( mesh.E1 );
            F2 = mesh.rotate_vf( F1 );
            
            I = repmat(1:mesh.nf,1,3); I = I(:); 
            J = (1:mesh.nf)'; J = [J; J+mesh.nf; J+2*mesh.nf];
            B1 = sparse(I,J,F1,mesh.nf,3*mesh.nf);
            B2 = sparse(I,J,F2,mesh.nf,3*mesh.nf);
            
            EB = [B1; B2];
            EBI = EB';
        end
    end
    
    methods (Static)
        
        function [ nf ] = normalize_f( f )
            nf = (f - min(f)) / (max(f)-min(f));
%             nf = f / (max(f)-min(f));
        end
        
        function [ snv ] = snormv( v )
            snv = sum(v.^2,2);
        end
        
        function [ nv ] = normv( v )
            nv = sqrt(sum(v.^2,2));
        end
        
        function [ nnv ] = normalize_vf( v )
            v = reshape(v,[],3);
            vn = MESH.normv(v); I = vn > 0;
            nnv = v;
            nnv(I,:) = v(I,:) ./ repmat(vn(I),1,3);
        end

        function [ B, BI, D ] = func_basis( M, k )
            L = - M.D * M.G;
            A = M.VA;
            W = A*L;
            tic; fprintf('\tEIGS of Laplace--Beltrami operator: ');
            [ev,el] = eigs(W,A,k,'SM');
            fprintf('%f sec\n',toc);
            
            [~,ii] = sort(diag(el)); el = el(ii,ii); ev = ev(:,ii);
            B = ev; BI = ev'*A; D = el;
        end
        
    end
end