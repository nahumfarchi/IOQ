function [thetas, periods, E, R] = solve_MIQ(A, b, thetas, periods, theta_tags, period_tags, direct_rounding)
    % function [thetas, periods, E] = solve_MIQ(A, b, thetas, periods, theta_tags, period_tags, direct_rounding)
    % 
    % Solve the given MIQ system by the normal equation and direct
    % rounding.
    %
    % Example:
    %   [local_frames, r] = create_local_frames(mesh);
    %   [A, b, t, p, t_tags, p_tags] = create_MIQ_system(m, N, f0, v0, local_frames, r); 
    %   [t, p, E] = solve_MIQ(A, b, t, p, t_tags, p_tags);
    
    if nargin < 7, direct_rounding = true; end
   
    n_thetas = sum(theta_tags>0);
    n_periods = sum(period_tags>0);
    
    I0 = diag([zeros(n_thetas, 1); ones(n_periods, 1)]);
    x = (I0+A'*A) \ (A'*b);
    %x = (A'*A) \ (A'*b);

    for fid = 1:length(thetas)
        if theta_tags(fid) > 0
            thetas(fid) = x(theta_tags(fid));
        end
    end

    for eid = 1:length(periods)
        if period_tags(eid) > 0
            if direct_rounding
                x(period_tags(eid)) = round(x(period_tags(eid)));
            end
            periods(eid) = x(period_tags(eid));
        end
    end

    E = norm(A*x-b)^2;

end

