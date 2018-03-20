function [exitflag] = check_grad(f, fcon, n, EPS)
    %CHECK_GRAD Check numerically the gradient of the given function.
    %
    % Input:
    %   TODO
    %   n - problem dimension

    if nargin < 4
        EPS = 1e-7;
    end

    %f = @(x) E(x, A, ba, B, bb);
    %fcon = @(x) deal([], C(:,1:end-1)*x-C(:,end));
    options = optimoptions(@fmincon, 'Algorithm', 'interior-point',...
    'CheckGradients',true, 'SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient',false,'MaxFunctionEvaluations',1);

    [~, ~, exitflag] = fmincon(f, zeros(n,1), [], [], [], [], [], [], fcon, options);
    
    %x = rand(n, 1);
    %[y,g] = f(x);
    %
    %for i = 1:n
    %%i = floor(n*rand)+1;
    %    x1 = x - EPS*delta(i, n);
    %    x2 = x + EPS*delta(i, n);
    %    numg = (f(x2) - f(x1)) / (2*EPS);
    %    assert(check_norm('g(i)', 'numg', 1e-7));
    %end

end

