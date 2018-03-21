classdef COREMESH < handle
    
    properties (Access='public')
        p       % parameters
        
        name1,name2
        M1,M2
        M1Q,M2Q
        
        TA1,TA2
        M1BIIf2v,M2BIIf2v
        
        M1B,M1BI,M1D    % Laplace--Beltrami eigen vectors and values on M1
        M2B,M2BI,M2D    % Laplace--Beltrami eigen vectors and values on M2
        
        map12,C12,k1    % map and fmap from M1 to M2
        map21,C21,k2    % map and fmap from M2 to M1
        usebasis
        
        GF1,GFT1
        GF2,GFT2
        m
        DOTSUM1,DOTSUM2
        
        HS,HL
        
        ew1,ew2         % curvature directions on M1 and M2
        WD1,WD2         % curvature weights on M1 and M2
        WDTA1,WDTA2
        cdv1,cdv2
        so1,so2
        
        n               % n-direction fields, currently n=4
        lp,ll           % weights for power and alignment energies
        maxdim          % max. dimension for solving with a direct method
        opt,opttool     % optimization tool used to minimize Eb and Ep
        abseps
        
        xf1,pf1         % smooth cross and power vector field on M1
        xf2,pf2         % smooth cross and power vector field on M2
        uv1,fuv1        % MIQ parametrization on M1
        uv2,fuv2        % MIQ parametrization on M2
        gradsz1,gradsz2
        
        E,ES,EC,EL      % energy values during minimization
        mtime,miter     % computation time and iteration count
        
        miqdir,miqbin
        qexdir,qexbin
        
        datadir
        outdir
        shareddir
        objsuffix
    end
    
    methods
        % Constructor :)
        function [ crm ] = COREMESH( name1, name2 )
            
            % parameters
            crm.p = inputParser;
            
            % fmaps & fvfs
            addOptional(crm.p,'map12',[]);
            addOptional(crm.p,'map21',[]);
            addOptional(crm.p,'k1',50);
            addOptional(crm.p,'k2',51);
            addOptional(crm.p,'usebasis',1);
            
            % cross fields
            addOptional(crm.p,'n',4);
            addOptional(crm.p,'cdv1',[]);
            addOptional(crm.p,'cdv2',[]);
            addOptional(crm.p,'so1',.1);
            addOptional(crm.p,'so2',.1);
            
            % optimization parameters
            addOptional(crm.p,'lp',.1);
            addOptional(crm.p,'ll',.1);
            addOptional(crm.p,'maxdim',2e4);
            addOptional(crm.p,'opt','direct');
            addOptional(crm.p,'opttool','minfunc');
            addOptional(crm.p,'abseps',1e-12);

            % libigl & libQEx
            addOptional(crm.p,'miqdir','../external/libigl/bin/build/');
            addOptional(crm.p,'miqbin','MIQ_bin');
            addOptional(crm.p,'qexdir','../external/libQEx/bin/');
            addOptional(crm.p,'qexbin','QEX_bin');
            addOptional(crm.p,'gradsz1',75);
            addOptional(crm.p,'gradsz2',75);
            
            % data directories
            addOptional(crm.p,'datadir','../data/');
            addOptional(crm.p,'outdir','experiments/');
            addOptional(crm.p,'shareddir','../');
            addOptional(crm.p,'objsuffix','');
            
            % create MESH data structures M1 and M2
            fprintf('Construct data structures\n\n');
            crm.name1 = name1; crm.M1 = MESH( crm.name1 ); 
            crm.name2 = name2; crm.M2 = MESH( crm.name2 );
        end
        
        % Problem parameters :)
        function set_param( crm, varargin )
            
            % parse & store input parameters
            parse( crm.p, varargin{:} );
            
            % fmaps & fvfs
            crm.map12 = crm.p.Results.map12;
            crm.map21 = crm.p.Results.map21;
            crm.k1 = crm.p.Results.k1;
            crm.k2 = crm.p.Results.k2;
            crm.usebasis = crm.p.Results.usebasis;
            
            % cross fields
            crm.n = crm.p.Results.n;
            crm.cdv1 = crm.p.Results.cdv1;
            crm.cdv2 = crm.p.Results.cdv2;
            crm.so1 = crm.p.Results.so1;
            crm.so2 = crm.p.Results.so2;
            
            % optimization parameters
            crm.lp = crm.p.Results.lp;
            crm.ll = crm.p.Results.ll;
            crm.abseps = crm.p.Results.abseps;
            crm.maxdim = crm.p.Results.maxdim;
            crm.opt = crm.p.Results.opt;
            crm.opttool = crm.p.Results.opttool;
            
            % libigl & libQEx
            crm.miqdir = crm.p.Results.miqdir;
            if ispc == 1
                crm.miqdir = '../external/libigl/bin/build_win64/';
            end
            crm.miqbin = crm.p.Results.miqbin;
            crm.gradsz1 = crm.p.Results.gradsz1;
            crm.gradsz2 = crm.p.Results.gradsz2;
            
            crm.qexdir = crm.p.Results.qexdir;
            if ispc == 1
                crm.qexdir = '../external/libQEx/bin_win64/';
            end
            crm.qexbin = crm.p.Results.qexbin;
            
            % data directories
            crm.datadir = crm.p.Results.datadir;
            crm.outdir = crm.p.Results.outdir;
            crm.shareddir = crm.p.Results.shareddir;
            crm.objsuffix = crm.p.Results.objsuffix;
            
            % 0. set preliminaries: meshes, eigenbases and fmaps
            crm.prep_prelim();
            fprintf('\n');
        end
        
        % Pre-process several computations :)
        function prep_prelim( crm )
            CM1 = crm.M1; 
            CM2 = crm.M2;
            
            crm.TA1 = spdiags([CM1.ta;CM1.ta],0,2*CM1.nf,2*CM1.nf);
            crm.TA2 = spdiags([CM2.ta;CM2.ta],0,2*CM2.nf,2*CM2.nf);
            
            % Laplace--Beltrami eigenfunctions and eigenvalues
            [crm.M1B,crm.M1BI,crm.M1D] = COREMESH.func_basis(CM1,crm.k1);
            [crm.M2B,crm.M2BI,crm.M2D] = COREMESH.func_basis(CM2,crm.k2);
            
            % mappings and functional mappings
            crm.prep_maps();
            
            % E_b and E_p useful matrices
            if crm.usebasis == 1
                crm.M1BIIf2v = crm.M1BI * crm.M1.If2v;
                crm.M2BIIf2v = crm.M2BI * crm.M2.If2v;
            else
                crm.M1BIIf2v = crm.M1.If2v;
                crm.M2BIIf2v = crm.M2.If2v;
            end
            
            % test functions
            f1 = crm.M1B*crm.M1D; f1 = f1(:,2:end);
            f2 = crm.M2B*crm.M2D; f2 = f2(:,2:end);
            if crm.usebasis == 1
                f1t = crm.M2B*(crm.C12*crm.M1D); f1t = f1t(:,2:end);
                f2t = crm.M1B*(crm.C21*crm.M2D); f2t = f2t(:,2:end);
            else
                f1t = crm.C12*f1;
                f2t = crm.C21*f2;
            end
            
            crm.m = min(crm.k1-1,crm.k2-1);
            
            w1 = diag(crm.M1D); w1 = [1; 1./w1(2:end)];
            w2 = diag(crm.M2D); w2 = [1; 1./w2(2:end)];
            
            crm.GF1 = CM1.G*f1;     crm.GF2 = CM2.G*f2;
            crm.GFT1 = CM2.G*f1t;   crm.GFT2 = CM1.G*f2t;
            for i = 1:crm.m
                crm.GF1(:,i) = w1(i)*COREMESH.get_pfld(CM1,crm.GF1(:,i),CM1.F1,crm.n);
                crm.GFT1(:,i) = w2(i)*COREMESH.get_pfld(CM2,crm.GFT1(:,i),CM2.F1,crm.n);
                
                crm.GF2(:,i) = w2(i)*COREMESH.get_pfld(CM2,crm.GF2(:,i),CM2.F1,crm.n);
                crm.GFT2(:,i) = w1(i)*COREMESH.get_pfld(CM1,crm.GFT2(:,i),CM1.F1,crm.n);
            end
            crm.GF1 = crm.GF1(:,1:crm.m); crm.GFT1 = crm.GFT1(:,1:crm.m);
            crm.GF2 = crm.GF2(:,1:crm.m); crm.GFT2 = crm.GFT2(:,1:crm.m);
            
            % aux. matrices
            II1 = repmat((1:CM1.nf)',3,1);
            JJ1 = 1:3*CM1.nf;
            SS1 = ones(numel(JJ1(:)),1);
            crm.DOTSUM1 = sparse(II1,JJ1,SS1,CM1.nf,3*CM1.nf);
            
            II2 = repmat((1:CM2.nf)',3,1);
            JJ2 = 1:3*CM2.nf;
            SS2 = ones(numel(JJ2(:)),1);
            crm.DOTSUM2 = sparse(II2,JJ2,SS2,CM2.nf,3*CM2.nf);
            
            crm.E = []; crm.ES = []; crm.EC = []; crm.EL = []; 
        end
        
        % Construct functional maps from M1 to M2 and vice versa :)
        function prep_maps( crm )
            CM1 = crm.M1; 
            CM2 = crm.M2;
            
            if isempty( crm.map21 ) == 1        % 1-1 correspondence
                crm.map12 = (1:CM1.nv)';
                crm.map21 = (1:CM2.nv)';
                
                crm.C12 = spdiags(ones(CM1.nv,1),0,CM1.nv,CM1.nv);
                crm.C21 = spdiags(ones(CM2.nv,1),0,CM2.nv,CM2.nv);
            end
            
            if size( crm.map21, 1 ) == CM1.nv
                if size( crm.map12, 2 ) == 4        % precise maps
                    
                    T2M = CM2.triangles(crm.map21(:,1),:);
                    W2M = crm.map21(:,2:4);
                    
                    T1M = CM1.triangles(crm.map12(:,1),:);
                    W1M = crm.map12(:,2:4);
                    
                    if crm.usebasis == 1
                        B2 = repmat(W2M(:,1),1,crm.k2) .* crm.M2B(T2M(:,1),:) + ...
                             repmat(W2M(:,2),1,crm.k2) .* crm.M2B(T2M(:,2),:) + ...
                             repmat(W2M(:,3),1,crm.k2) .* crm.M2B(T2M(:,3),:) ;
                        crm.C12 = B2 \ crm.M1B;
                    
                        B1 = repmat(W1M(:,1),1,crm.k1) .* crm.M1B(T1M(:,1),:) + ...
                             repmat(W1M(:,2),1,crm.k1) .* crm.M1B(T1M(:,2),:) + ...
                             repmat(W1M(:,3),1,crm.k1) .* crm.M1B(T1M(:,3),:) ;
                        crm.C21 = B1 \ crm.M2B;
                    else
                        II1 = repmat((1:CM2.nv)',3,1);
                        crm.C12 = sparse(II1,double(T1M(:)),W1M(:),CM2.nv,CM1.nv);
                        
                        II2 = repmat((1:CM1.nv)',3,1);
                        crm.C21 = sparse(II2,double(T2M(:)),W2M(:),CM1.nv,CM2.nv);
                    end
                else                                  % vertex to vertex maps
                    if crm.usebasis == 1
                        crm.C12 = crm.M2B(crm.map21,:) \ crm.M1B;
                        crm.C21 = crm.M1B(crm.map12,:) \ crm.M2B;
                    else
                        crm.C12 = sparse(1:CM2.nv,crm.map12,ones(CM2.nv,1),CM2.nv,CM1.nv);
                        crm.C21 = sparse(1:CM1.nv,crm.map21,ones(CM1.nv,1),CM1.nv,CM2.nv);
                    end
                end
            else                                      % sparse correspondences
%                 [crm.C12,crm.C21] = compute_fmap([0;1./diag(crm.M1D(2:end,2:end))], crm.M1B, sum(CM1.ta), crm.map12, [0;1./diag(crm.M2D(2:end,2:end))], crm.M2B, sum(CM2.ta), crm.map21);
                [crm.map12,crm.C12] = compute_fmap2(CM1,crm.M1B,crm.M1D,crm.k1,CM2,crm.M2B,crm.M2D,crm.k2,crm.map12);
                [crm.map21,crm.C21] = compute_fmap2(CM2,crm.M2B,crm.M2D,crm.k2,CM1,crm.M1B,crm.M1D,crm.k1,crm.map21);
                
                crm.map12 = crm.map12-1; crm.map21 = crm.map21-1;
                dlmwrite([crm.shareddir crm.name1 '_to_' crm.name2 '.fmap.map'],crm.map12,' ');
                dlmwrite([crm.shareddir crm.name2 '_to_' crm.name1 '.fmap.map'],crm.map21,' ');
                crm.map12 = crm.map12+1; crm.map21 = crm.map21+1;
            end
%             crm.dbg_maps();
        end
        
        function [ xfm ] = push_xfld_p2p( crm, xf )
            
            D = COREMESH.map_diff(crm.M1,crm.M2);
                        
            xfm = zeros(crm.M2.nf,3);
            for i = 1:crm.M2.nf
                d = reshape( D(i,:), 3, 3 );
                xfm(i,:) = (d*xf(i,:)')';
            end
            xfm = MESH.normalize_vf( xfm );
        end
        
        % Curvature directions and weights :)
        function comp_shapeop( crm )
            CM1 = crm.M1; 
            CM2 = crm.M2;
            
            tic
            if ~crm.ll == 0
                if isempty(crm.cdv1) == 1
                    [w1,wd1] = SHAPEOP.shape_operator( CM1 );
                    fprintf('name=%s, Max. curvature b4 clamp: %g\n', ...
                            CM1.name, max(wd1));
                    wd1(wd1 < crm.so1) = 0;
                else
                    w1 = crm.cdv1; 
                    wd1 = MESH.normv( reshape(w1,[],3) );
                end
                if isempty(crm.cdv2) == 1
                    [w2,wd2] = SHAPEOP.shape_operator( CM2 );
                    fprintf('name=%s, Max. curvature b4 clamp: %g\n', ...
                            CM2.name, max(wd2));
                    wd2(wd2 < crm.so2) = 0;
                else
                    w2 = crm.cdv2;
                    wd2 = MESH.normv( reshape(w2,[],3) );
                end
                
                w1 = COREMESH.get_pfld( CM1, w1, CM1.F1, crm.n );
                w2 = COREMESH.get_pfld( CM2, w2, CM2.F1, crm.n );
                
                crm.ew1 = CM1.EB * w1(:);
                crm.ew2 = CM2.EB * w2(:);
                
                crm.WD1 = spdiags(repmat(wd1,2,1),0,2*CM1.nf,2*CM1.nf);
                crm.WD2 = spdiags(repmat(wd2,2,1),0,2*CM2.nf,2*CM2.nf);
            else
                crm.ew1 = zeros(2*CM1.nf,1);
                crm.ew2 = zeros(2*CM2.nf,1);
                crm.WD1 = 0; crm.WD2 = 0;
            end
            
            crm.WDTA1 = crm.WD1*crm.TA1;
            crm.WDTA2 = crm.WD2*crm.TA2;
            fprintf('\tComputing power version and weights: %f sec\n',toc);
%             crm.dbg_shapeop();
        end
        
        % Align the smooth power vector fields on M1 and M2
        function [ e, g ] = comp_energy( crm, ey, ~ )
            CM1 = crm.M1; Bfv1 = crm.M1BIIf2v; 
            CM2 = crm.M2; Bfv2 = crm.M2BIIf2v;
            
            ey1 = ey(1:2*CM1.nf); ey2 = ey(2*CM1.nf+1:end);
            
            % smoothness terms
            Ley1 = CM1.VLn*ey1;
            Ley2 = CM2.VLn*ey2;
            es = .5*ey1'*Ley1 + .5*ey2'*Ley2;
            
            % consistency terms
            ry1 = repmat(CM1.EBI*ey1,1,crm.m);
            ry2 = repmat(CM2.EBI*ey2,1,crm.m);
            
            YGF1 = crm.DOTSUM1*(ry1.*crm.GF1);
            YGFT1 = crm.DOTSUM2*(ry2.*crm.GFT1);
            
            YGF2 = crm.DOTSUM2*(ry2.*crm.GF2);
            YGFT2 = crm.DOTSUM1*(ry1.*crm.GFT2);
            
            ct1 = Bfv1*YGF1 - crm.C21*(Bfv2*YGFT1);
            ct2 = Bfv2*YGF2 - crm.C12*(Bfv1*YGFT2);
            ec = .5*sum( dot(ct1,ct1,1) ) + .5*sum( dot(ct2,ct2,1) );
            
            % alignment terms
            yw1 = ey1-crm.ew1; WDTAyw1 = crm.WDTA1*yw1; 
            yw2 = ey2-crm.ew2; WDTAyw2 = crm.WDTA2*yw2; 
            el = .5*yw1'*WDTAyw1 + .5*yw2'*WDTAyw2;
            
            % energy
            a = crm.lp; b = crm.ll;
            e = (1-b)*((1-a)*es + a*ec) + b*el;
            
            if nargout > 1
                % smoothness terms
                gs = [Ley1; Ley2];
                
                % consistency term
                bct1 = crm.GF1 .* repmat(Bfv1'*ct1,3,1);
                bctt1 = crm.GFT1 .* repmat(Bfv2'*(crm.C21'*ct1),3,1);
                
                bct2 = crm.GF2 .* repmat(Bfv2'*ct2,3,1);
                bctt2 = crm.GFT2 .* repmat(Bfv1'*(crm.C12'*ct2),3,1);
                
                gc = [ sum( CM1.EBI'*(bct1 - bctt2), 2 ); ...
                       sum( CM2.EBI'*(bct2 - bctt1), 2 ) ];
                
                % alignment terms
                gl = [WDTAyw1; WDTAyw2];
                
                % gradient
                g = (1-b)*((1-a)*gs + a*gc) + b*gl;
                
%                 % DEBUG
%                 edbg = .5*ey'*(1-b)*((1-a)*gs+a*gc) + .5*b*[yw1; yw2]'*gl;
%                 fprintf('e=%g, edbg=%g, |e-edbg|=%g\n',e,edbg,abs(e-edbg));
            end
            
            % stats
            crm.ES = [crm.ES; es]; crm.EC = [crm.EC; ec]; crm.EL = [crm.EL; el];
        end
        
        function [ g ] = comp_gradient( crm, ey, ~ )
            [~,g] = crm.comp_energy( ey );
        end
        
        function [ Hv ] = comp_hessianv( crm, ev, ey, ~ )
            CM1 = crm.M1; Bfv1 = crm.M1BIIf2v; 
            CM2 = crm.M2; Bfv2 = crm.M2BIIf2v;
            
            a = crm.lp; b = crm.ll;
            Hv = ((1-b)*(1-a)*crm.HS + b*crm.HL)*ev;
            
            ev1 = ev(1:2*CM1.nf); 
            ev2 = ev(2*CM1.nf+1:end);
            
            % consistency terms
            ry1 = repmat(CM1.EBI*ev1,1,crm.m);
            ry2 = repmat(CM2.EBI*ev2,1,crm.m);
            
            YGF1 = crm.DOTSUM1*(ry1.*crm.GF1);
            YGFT1 = crm.DOTSUM2*(ry2.*crm.GFT1);
            
            YGF2 = crm.DOTSUM2*(ry2.*crm.GF2);
            YGFT2 = crm.DOTSUM1*(ry1.*crm.GFT2);
            
            ct1 = Bfv1*YGF1 - crm.C21*(Bfv2*YGFT1);
            ct2 = Bfv2*YGF2 - crm.C12*(Bfv1*YGFT2);
            
            bct1 = crm.GF1 .* repmat(Bfv1'*ct1,3,1);
            bctt1 = crm.GFT1 .* repmat(Bfv2'*(crm.C21'*ct1),3,1);

            bct2 = crm.GF2 .* repmat(Bfv2'*ct2,3,1);
            bctt2 = crm.GFT2 .* repmat(Bfv1'*(crm.C12'*ct2),3,1);
            gc = [ sum( CM1.EBI'*(bct1 - bctt2), 2 ); ...
                   sum( CM2.EBI'*(bct2 - bctt1), 2 ) ];
            
            Hv = Hv + (1-b)*a*gc;
        end
        
        function [ Hv ] = comp_hessianv_qp( crm, B, ev )
            Hv = zeros( size(ev) );
            % consistency term
            for i = 1:size(ev,2)
                Hv(:,i) = crm.comp_hessianv( ev(:,i), [] );
            end
        end
        
        function [ B, f, l, u ] = quadprogParams( crm )
            CM1 = crm.M1; CM2 = crm.M2; 
            a = crm.lp; b = crm.ll;
            B = (1-b)*(1-a)*crm.HS + b*crm.HL;
            f = - b*[crm.ew1'*crm.WDTA1 crm.ew2'*crm.WDTA2]';
            l = -1000*ones(2*CM1.nf+2*CM2.nf,1);
            u = 1000*ones(2*CM1.nf+2*CM2.nf,1);
        end
        
        function opt_power_vfs( crm )
            CM1 = crm.M1; CM2 = crm.M2;
            ze1 = spalloc(2*CM1.nf,2*CM2.nf,1);
            ze2 = spalloc(2*CM2.nf,2*CM1.nf,1);
            crm.HS = [CM1.VLn ze1; ze2 CM2.VLn]; 
            crm.HL = [crm.WDTA1 ze1; ze2 crm.WDTA2];
            
            tic;
            [ey1,~] = eigs(CM1.VLn,[],6,'SM'); ey1 = ey1(:,1);
            [ey2,~] = eigs(CM2.VLn,[],6,'SM'); ey2 = ey2(:,1);
            ey = [ey1; ey2];
            
            if strcmp(crm.opt,'direct') == 1
                
                opts.directparams = @crm.quadprogParams;
                opts.TolFun = crm.abseps;
                CH = @(B,y) crm.comp_hessianv_qp( B, y );
                [epf,output] = OPT.direct(CH,ey,opts);
                
            else
                opts = [];
                opts.display = 'iter';
                opts.optTol = crm.abseps; 
                opts.progTol = crm.abseps;
                opts.MaxFunEvals = 1000;
                opts.MaxIter = 1000;
%                 opts.Method = 'pnewton0';
%                 opts.HvFunc = @(g,y,varargin) crm.comp_hessianv( g, y, varargin );
                [epf,output] = OPT.iter(@crm.comp_energy,...
                                    @crm.comp_gradient,...
                                    crm.opttool,ey,opts);
                crm.E = output.trace.fval;
            end
            crm.mtime = toc;
            fprintf('\tComputation took: %f sec\n',crm.mtime);
            
            % stats    
            crm.miter = output.iterations;
            
            % extract y_1 and y_2 and the corresponding x_1 and x_2
            crm.pf1 = reshape(CM1.EBI*epf(1:2*CM1.nf),[],3);
            crm.pf2 = reshape(CM2.EBI*epf(2*CM1.nf+1:end),[],3);
            crm.xf1 = COREMESH.get_xfld( CM1, crm.pf1, CM1.F1, crm.n );
            crm.xf2 = COREMESH.get_xfld( CM2, crm.pf2, CM2.F1, crm.n );
            crm.xf1 = reshape(crm.xf1,[],3);
            crm.xf2 = reshape(crm.xf2,[],3);
            
            % store results in files
            fid = fopen( [ crm.outdir CM1.name '_nrosy.txt' ], 'w' );
            fprintf(fid, '%g %g %g\n',crm.xf1'); fclose(fid);
            
            fid = fopen( [ crm.outdir CM2.name '_nrosy.txt' ], 'w' );
            fprintf(fid, '%g %g %g\n',crm.xf2'); fclose(fid);
        end
        
        % Execute libigl to compute a parametrization using comiso :)
        function [ uv, fuv ] = comp_uv( crm, mesh, gsz )
            
            % output goes to outdir with names: meshname '_uv' or '_fuv'
            eval( sprintf('!"%s%s" %s %s %s %d',crm.miqdir,crm.miqbin,...
                  crm.datadir,mesh.name,crm.outdir,gsz ) );
            
            % load UV and FUV
            uv = load([crm.outdir mesh.name '_uv.txt']);         % size 2d_nv x 2
            fuv = load([crm.outdir mesh.name '_fuv.txt']) + 1;   % size nf x 3
        end
        
        % Store the resulting texture coordinates in OBJ :)
        function store_obj( crm, mesh, uv, fuv )
            
            options.nm_file = [crm.shareddir 'cross.png'];
            options.face_texcorrd = fuv;
            options.object_texture = uv;
            write_obj( crm.outdir, [mesh.name crm.objsuffix '.obj'], ...
                       mesh.vertices, mesh.triangles, options );
                   
        end
        
        % Extract quadrangular mesh from the parametrization :)
        function [ MQ ] = comp_quad( crm, name )
            triname = [ crm.outdir name crm.objsuffix '.obj' ];
            quadname = [ crm.outdir name '_quad' crm.objsuffix '.obj' ];
            
            runqex = sprintf('"%s%s" %s %s',...
                             crm.qexdir,crm.qexbin,triname,quadname);
            if ispc == 1
                eval( ['!' runqex] );
            else
                chgld = ['!export DYLD_LIBRARY_PATH=' crm.qexdir];
                eval( [chgld '; ' runqex] );
            end
            
            % store into MQ
            [vq,fq] = read_obj(quadname);
            MQ.vertices = vq'; MQ.quads = fq';
        end
        
        % Main algorithm, including all the required sub-steps :)
        function solve( crm )
            
            % 1 compute curvature direction
            fprintf('Compute curvature directions\n\n');
            crm.comp_shapeop();
            fprintf('\n');

            % 2. optimize for smooth & consistent power vector fields
            fprintf('Optimize for smooth & consistent power vector fields\n\n');
            crm.opt_power_vfs();
            fprintf('\n');
            
            % 3. libigl: compute parametrizations for M1/2
            fprintf('libigl: compute parametrizations\n\n');
            [crm.uv1,crm.fuv1] = crm.comp_uv( crm.M1, crm.gradsz1 );
            [crm.uv2,crm.fuv2] = crm.comp_uv( crm.M2, crm.gradsz2 );
            
            % save results to OBJ files
            fprintf('Store results in OBJ files\n\n');
            crm.store_obj( crm.M1, crm.uv1, crm.fuv1 );
            crm.store_obj( crm.M2, crm.uv2, crm.fuv2 );
            
            % 4. libQEx: generate quad meshes
            fprintf('Generate quad meshes\n\n');
            crm.M1Q = crm.comp_quad( crm.M1.name );
            crm.M2Q = crm.comp_quad( crm.M2.name );
            
            % 5. plot the result
            figure; MESH_VIS.mesh( crm.M1Q );
            figure; MESH_VIS.mesh( crm.M2Q );
        end
        
        function gen_quad( crm, M, xf, gsz )
            fid = fopen( [ crm.outdir M.name '_nrosy.txt' ], 'w' );
            fprintf(fid, '%g %g %g\n',xf'); fclose(fid);
            
            fprintf('libigl: compute parametrization\n\n');
            [uv,fuv] = crm.comp_uv( M, gsz );
            
            fprintf('Store results in an OBJ file\n\n');
            crm.store_obj( M, uv, fuv );
            
            fprintf('Generate quad meshes\n\n');
            crm.comp_quad( M.name );
        end
        
        function [bf,pf] = design_xfld( crm, mesh )
            
            % mesh.F1/F2 are used for the frame
            [pf,bf] = godf_dirichlet(mesh,crm.n);
%             figure; MESH_VIS.vf(mesh,sf,'nRosy',crm.n);
            
            % store cross field in outdir
            fid = fopen( [ crm.outdir mesh.name '_nrosy.txt' ], 'w' );
            fprintf(fid, '%g %g %g\n',bf'); fclose(fid);
        end
        
        function [xfm,pfm] = push_xfld( crm, xf, pf )
            
            % using the differential of a p2p map
%             sfm = crm.push_xfld_p2p( xf );
            
            % using the functional differential of a p2p map
%             sfm = crm.push_xfld_fmap_angle( rf );
            xfm = crm.push_xfld_fmap_proj( pf );
            pfm = xfm;
            
            % store cross field in outdir
            fid = fopen( [ crm.outdir crm.M2.name '_nrosy.txt' ], 'w' );
            fprintf(fid, '%g %g %g\n',xfm'); fclose(fid);
        end
        
        function [ xfm ] = push_xfld_fmap_proj( crm, pf )
            % push M1.F1
            M2F1 = crm.push_xfld_p2p( crm.M1.F1 );
            M2F2 = reshape( crm.M2.R * M2F1(:), [], 3 );
            
            % push projections of rf onto the frame
            pu1 = dot( pf, crm.M1.F1, 2 ); pu1v = crm.M1.If2v*pu1;
            pv1 = dot( pf, crm.M1.F2, 2 ); pv1v = crm.M1.If2v*pv1;
            
            pu2v = crm.M2B * (crm.C12 * (crm.M1BI * pu1v));
            pv2v = crm.M2B * (crm.C12 * (crm.M1BI * pv1v));
            
            pu2 = crm.M2.Iv2f*pu2v;
            pv2 = crm.M2.Iv2f*pv2v;
                        
            % rotate RF2 in theta2 angle
            xfm = M2F1 .* repmat( pu2, 1, 3 ) + ...
                  M2F2 .* repmat( pv2, 1, 3 );
        end
        
        function [ xfm ] = push_xfld_fmap_angle( crm, pf )
            % push M1.F1
            M2F1 = crm.push_xfld_p2p( crm.M1.F1 );
            M2F2 = reshape( crm.M2.R * M2F1(:), [], 3 );
            
            % push angle
            a1 = atan2(dot(pf,crm.M1.F2,2),dot(pf,crm.M1.F1,2));
            
            a1v = crm.M1.If2v*a1;
            a2v = crm.M2B * (crm.C12 * (crm.M1BI * a1v));
            a2 = crm.M2.Iv2f*a2v;
%             a2 = a1;
                        
            % rotate RF2 in theta2 angle
            xfm = M2F1 .* repmat(cos(a2/crm.n),1,3)  + ...
                  M2F2 .* repmat(sin(a2/crm.n),1,3);
        end
        
        function solve_naive( crm )
            
            % set preliminaries
            fprintf('Construct data structures\n\n');
            crm.prep_prelim();
            fprintf('\n');
            
            % design & store a cross field on M1
            fprintf('Design cross field on %s\n\n', crm.name1);
            [crm.xf1,crm.pf1] = crm.design_xfld( crm.M1 );
            
            % push cross field to M2 & store
            fprintf('Push cross field to %s\n\n', crm.name2);
            [crm.xf2,crm.pf2] = crm.push_xfld( crm.xf1, crm.pf1 );
            
            % libigl: compute parametrizations for M1/2
            fprintf('libigl: compute parametrizations\n\n');
            [crm.uv1,crm.fuv1] = crm.comp_uv( crm.M1 );
            [crm.uv2,crm.fuv2] = crm.comp_uv( crm.M2 );
            
            % save results to OBJ files
            fprintf('Store results in OBJ files\n\n');
            crm.store_obj( crm.M1, crm.uv1, crm.fuv1 );
            crm.store_obj( crm.M2, crm.uv2, crm.fuv2 );
        end
        
        function dbg_maps( crm )
            % map from M1 to M2, mapi from M2 to M1
            f1 = crm.M1B(:,10);
            f2 = crm.M2B(:,20);

            if size( crm.map12, 1 ) == crm.M2.nv && ...
               size( crm.map21, 1 ) == crm.M1.nv
               
                if size( crm.map12, 2 ) == 4
                    T1M = crm.M1.triangles(crm.map12(:,1),:); 
                    T2M = crm.M2.triangles(crm.map21(:,1),:); 
                    W1M = crm.map12(:,2:4);
                    W2M = crm.map21(:,2:4);

                    f12 = W1M(:,1) .* f1(T1M(:,1),:) + ...
                          W1M(:,2) .* f1(T1M(:,2),:) + ...
                          W1M(:,3) .* f1(T1M(:,3),:) ;

                    f21 = W2M(:,1) .* f2(T2M(:,1),:) + ...
                          W2M(:,2) .* f2(T2M(:,2),:) + ...
                          W2M(:,3) .* f2(T2M(:,3),:) ;
                else
                    f12 = f1(crm.map12);
                    f21 = f2(crm.map21);
                end
                
                figure; MESH_VIS.func(crm.M1,f1);
                figure; MESH_VIS.func(crm.M2,f12);
                figure; MESH_VIS.func(crm.M2,f2);
                figure; MESH_VIS.func(crm.M1,f21);
            end
            
            % C12 from L2(M1) to L2(M2), C21 from L2(M2) to L2(M1)
            if size(crm.C12,1) == crm.M2.nv
                g1 = f1; g12 = crm.C12*g1;
                g2 = f2; g21 = crm.C21*g2;
            else
                g1 = crm.M1BI*f1; g12 = crm.M2B*crm.C12*g1; g1 = crm.M1B*g1;
                g2 = crm.M2BI*f2; g21 = crm.M1B*crm.C21*g2; g2 = crm.M2B*g2;
            end
            
            figure; MESH_VIS.func(crm.M1,g1);
            figure; MESH_VIS.func(crm.M2,g12);
            figure; MESH_VIS.func(crm.M2,g2);
            figure; MESH_VIS.func(crm.M1,g21);
        end
        
        function dbg_shapeop( crm )
            xw1 = COREMESH.get_xfld(crm.M1,reshape(crm.M1.EBI*crm.ew1,[],3),crm.M1.F1,crm.n);
            xw2 = COREMESH.get_xfld(crm.M2,reshape(crm.M2.EBI*crm.ew2,[],3),crm.M2.F1,crm.n);
            figure; MESH_VIS.vf(crm.M1,xw1,'nRosy',2);
            figure; MESH_VIS.vf(crm.M2,xw2,'nRosy',2);
            figure; MESH_VIS.func(crm.M1,diag(crm.WD1(1:crm.M1.nf,1:crm.M1.nf)));
            figure; MESH_VIS.func(crm.M2,diag(crm.WD2(1:crm.M2.nf,1:crm.M2.nf)));
        end
        
        function plot_energy_terms( crm )
            a = crm.lp; b = crm.ll;
            
            figure; 
            plot((1-b)*(1-a)*crm.ES); hold on; 
            plot((1-b)*a*crm.EC); 
            plot(b*crm.EL);
            legend('Es','Ec','El');
            hold off;

        end
        
        function [ ec ] = eval_consistency( crm, x1, x2 )
            CM1 = crm.M1; Bfv1 = CM1.If2v;
            CM2 = crm.M2; Bfv2 = CM2.If2v;
            
            y1 = COREMESH.get_pfld(CM1,x1,CM1.F1,crm.n);
            y2 = COREMESH.get_pfld(CM2,x2,CM2.F1,crm.n);
            
            if size( crm.map12, 2 ) == 4
                    
                T2M = CM2.triangles(crm.map21(:,1),:);
                W2M = crm.map21(:,2:4);

                T1M = CM1.triangles(crm.map12(:,1),:);
                W1M = crm.map12(:,2:4);

                II1 = repmat((1:CM2.nv)',3,1);
                CT12 = sparse(II1,double(T1M(:)),W1M(:),CM2.nv,CM1.nv);

                II2 = repmat((1:CM1.nv)',3,1);
                CT21 = sparse(II2,double(T2M(:)),W2M(:),CM1.nv,CM2.nv);
            else                                  
                CT12 = sparse(1:CM2.nv,crm.map12,ones(CM2.nv,1),CM2.nv,CM1.nv);
                CT21 = sparse(1:CM1.nv,crm.map21,ones(CM1.nv,1),CM1.nv,CM2.nv);
            end
            
%             % choose szidx random vertices
%             szidx = 1000; idx = randi(CM.nv,szidx);
%             
%             ec = zeros(CM.nv,1); tic
%             for i = 1:szidx
%                 f = zeros(CM.nv,1); f(idx(i)) = 1;
%                 
%                 gf = COREMESH.get_pfld(CM,CM.G*f,CM.F1,4);
%                 gft= COREMESH.get_pfld(CM,CM.G*(CT*f),CM.F1,4);
%             
%                 ygf = DS*(y .* gf);
%                 ygft = DS*(y .* gft);
%                 diff = CM.If2v*ygf-CT*(CM.If2v*ygft);
%                 ec(idx(i)) = .5*(diff'*CM.VA*diff);
%             end
%             ec = sum(ec);
            ec = 0;
            
            % additionally, use the LB eigenfunctions
            ry1 = repmat(y1,1,crm.m);
            ry2 = repmat(y2,1,crm.m);
            
            YGF1 = crm.DOTSUM1*(ry1.*crm.GF1);
            YGFT1 = crm.DOTSUM2*(ry2.*crm.GFT1);
            
            YGF2 = crm.DOTSUM2*(ry2.*crm.GF2);
            YGFT2 = crm.DOTSUM1*(ry1.*crm.GFT2);
            
            ct1 = Bfv1*YGF1 - CT21*(Bfv2*YGFT1);
            ct2 = Bfv2*YGF2 - CT12*(Bfv1*YGFT2);
            
            for i = 1:crm.m
                ec = ec + .5*ct1(:,i)'*CM1.VA*ct1(:,i) ...
                        + .5*ct2(:,i)'*CM2.VA*ct2(:,i);
            end
            
            toc
        end
    end
    
    methods (Static)
        
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
        
        function [ D, MD ] = map_diff( M1, M2 )
            % assume same connectivity
            X1 = M1.vertices; T1 = M1.triangles; N1 = M1.Nf;
            X2 = M2.vertices; T2 = M2.triangles; N2 = M2.Nf;

            D = zeros(9,M2.nf);
            for i = 1:M2.nf
                
                x1 = [ X1(T1(i,:),:); X1(T1(i,1),:) + N1(i,:) ];
                ps1 = mean(x1);
                x1c = x1 - repmat(ps1,4,1);

                x2 = [ X2(T2(i,:),:); X2(T2(i,1),:) + N2(i,:) ];
                ps2 = mean(x2);
                x2c = x2 - repmat(ps2,4,1);
                
                d = x2c' * pinv(x1c');
                d = d'; D(:,i) = d(:);
            end
            
            II = repmat(1:3*M2.nf,3,1);
            JJ = reshape(1:3*M2.nf,[],3)'; JJ = repmat(JJ(:),3,1);
            SS = [ reshape(D(1:3,:),[],1); ...
                   reshape(D(4:6,:),[],1); ...
                   reshape(D(7:9,:),[],1) ];
            MD = sparse(II(:),JJ(:),SS(:),3*M2.nf,3*M2.nf);
        end
        
        function [ xf ] = get_xfld( M, pf, bf, n )
            rbf = reshape(M.R*bf(:),[],3);
            pf = MESH.normalize_vf( reshape(pf,[],3) );
            xf = pf;
            
            zi = MESH.normv( pf ) > eps;
            ti = atan2(dot(pf(zi,:),rbf(zi,:),2),...
                       dot(pf(zi,:),bf(zi,:),2)) / n;
            
            xf(zi,:) = repmat(cos(ti),1,3) .* bf(zi,:) + ...
                       repmat(sin(ti),1,3) .* rbf(zi,:);    
            xf = xf(:);
        end
        
        function [ pf ] = get_pfld( M, xf, bf, n )
            rbf = reshape(M.R*bf(:),[],3);
            xf = MESH.normalize_vf( reshape(xf,[],3) );
            pf = xf;
            
            zi = MESH.normv( xf ) > eps;
            ti = atan2( dot(xf(zi,:),rbf(zi,:),2),...
                        dot(xf(zi,:),bf(zi,:),2) ) * n;
            
            pf(zi,:) = repmat(cos(ti),1,3) .* bf(zi,:) + ...
                       repmat(sin(ti),1,3) .* rbf(zi,:);
            pf = pf(:);
        end
    end
end