function test_cycle_defects()

    fp = '../data/bunnies_small/bunny_res1.off';
    %fp = '../data/torus_s0.off';
    
    m = Mesh(fp); nv = m.nV;

    %% Standard TCODS cycle basis (i.e., a cycle around each vertex and the generators)
    vf_1ring = m.vf_1ring;
    tic; vert_cycles = m.vert_cycles; 
    vert_defects = m.cycle_defects(vert_cycles); toc
    for i = 1:nv
        grp1 = vf_1ring{i};
        grp2 = vert_cycles{i};
        assert(length(intersect(grp1, grp2)) == length(grp1))
        assert(grp2(1) == grp2(end))
    end

    tic; Kg = m.gaussian_cur; toc
    assert(norm(abs(Kg) - abs(vert_defects)) < 1e-10)

    gen_cycles = m.generator_cycles;
    cycles_tc = cell(numel(vert_cycles) + numel(gen_cycles), 1);
    cycles_tc(1:numel(vert_cycles)) = vert_cycles;
    cycles_tc(numel(vert_cycles)+1:end) = gen_cycles;
    Gamma_tc = m.dual_cycles_to_1forms(cycles_tc);
    
    %%
    T = pi / 2;
    [local_frames, r] = create_local_frames(m);
    EPS = 10.^(0:-1:-5);
    
    for i = 1:length(EPS)
        eps = EPS(i);
        lb = r - eps;
        ub = r + eps;
        POS = Gamma_tc; POS(POS<0) = 0;
        NEG = Gamma_tc; NEG(NEG>0) = 0;
        ac = ceil((POS'*lb - NEG'*ub) / T);
        bc = floor((POS'*ub - NEG'*lb) / T);
        cycle_width = bc - ac;
        
        disp(sum(cycle_width==0))
        
        cycle_widths{i} = cycle_width;
        acs{i} = ac;
        bcs{i} = bc;
        num_zeros(i) = sum(cycle_width==0);
        max_width(i) = max(cycle_width);
        min_width(i) = min(cycle_width);
    end
    figure; hold on
    plot(EPS, num_zeros, '--x') 
    plot(EPS, max_width, '--x')
    plot(EPS, min_width, '--x')
    legend('Number of zeros', 'max width', 'min width'); xlabel('eps')
    
    FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
    USE_GPU = false;

    % run IOQ
    rng(SEED)
    [alpha, beta, connection] = IOQ(...
        m.V, m.F, ...
        'GPU', USE_GPU);
    
    cmp=[cycle_widths{2}, acs{2}, bcs{2}, alpha];

    % run TCODS
    k = [alpha; beta];
    res_tc1 = TCODS(m, ...
        'k', k, ...
        'f0', FACE0, ...
        'theta0', THETA0, ...
        'degree', DEGREE, ...
        'CreateFField', true, ...
        'Duplicate', true, ...
        'gConstraintVec', GVEC, ...
        'constraints', [1, 0; 2, pi]);
    
    alpha2 = alpha; alpha2(9) = alpha2(9)-3; alpha2(16) = alpha2(16)+3;
    k2 = [alpha2; beta];
    res_tc2 = TCODS(m, ...
        'k', k2, ...
        'f0', FACE0, ...
        'theta0', THETA0, ...
        'degree', DEGREE, ...
        'CreateFField', true, ...
        'Duplicate', true, ...
        'gConstraintVec', GVEC);
    
    fprintf('E1 = %g, E2 = %g\n', res_tc1.miq_energy, res_tc2.miq_energy)
    

    %% Strictly fundamental (maybe?) dual basis cycles
    %[Gamma2, Ad2, cycles2] = SFCBasis(m2);
    [Gamma0, defects0, cycles0] = m.SFCBasisd();
    Gamma0_perm = reorder_cycle_basis(Gamma0); 
   
    assert(size(Gamma_tc, 1) == size(Gamma0, 1))
    assert(size(Gamma_tc, 2)-1 == size(Gamma0, 2))
    assert(norm(sortrows(Gamma0) - sortrows(Gamma0_perm), 'fro') < 1e-10)
    
    %% Strictly fundamental primal basis cycles
    [d0, d1] = get_exterior_derivatives(m);
    [Gamma1, defects1, cycles1] = m.SFCBasisp();
    Gamma1_perm = reorder_cycle_basis(Gamma1);
    
    figure;
    subplot(121); m.plotCycles(cycles0)
    subplot(122); m.plotCycles(cycles1)
    set(gcf, 'WindowStyle', 'docked')
    
    figure; 
    subplot(231); spy(d1'); title('d1^T')
    subplot(232); spy(Gamma1); title('Gamma1')
    subplot(233); spy(Gamma1_perm); title('Gamma1 perm')
    subplot(234); spy(Gamma_tc); title('Gamma tc')
    subplot(235); spy(Gamma0); title('Gamma0')
    subplot(236); spy(Gamma0_perm); title('Gamma0 perm')
    set(gcf, 'WindowStyle', 'docked')
    
    figure;
    subplot(131); spy(d1*Gamma_tc); title('d1*gamma tc')
    subplot(132); spy(d1*Gamma0); title('d1*Gamma0')
    subplot(133); spy(Gamma1'*Gamma0); title('Gamma1''*Gamma0')
    set(gcf, 'WindowStyle', 'docked')
    
    assert(norm(d1*Gamma0, 'fro') == 0)
    assert(norm(Gamma1'*Gamma0, 'fro') == 0);
    
function Gamma_perm = reorder_cycle_basis(Gamma)
    [ne, n_cycles] = size(Gamma);

    % Find edge that appear only in one cycle and move them to the start
    inds = find(sum(abs(Gamma), 2) == 1);
    ninds = setdiff(1:ne, inds);
    Gamma_perm = [Gamma(inds, :); Gamma(ninds, :)];

    % Reorder the first n_cycles rows so that the submatrix is equal to I
    [r,c] = find(Gamma_perm(1:n_cycles, :));
    % For each row, find the first non-zero column
    fi = accumarray(r, c, [n_cycles,1], @min, NaN);
    assert(~any(isnan(fi)))
    % Sort rows by earliest non-zero
    [~, si] = sort(fi);
    Gamma_perm = [Gamma_perm(si, :); Gamma_perm(n_cycles+1:end, :)];

    assert(norm(abs(Gamma_perm(1:n_cycles,:)) - speye(n_cycles), 'fro') < 1e-10)

    