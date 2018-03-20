function [pass, nrm_out] = check_norm(str1, str2, varargin)
%CHECK_NORM Summary of this function goes here
%   Detailed explanation goes here

    p = inputParser;
    addOptional(p, 'Tol', 1e-7);
    addOptional(p, 'Log', -1);
    parse(p, varargin{:});
    
    TOL = p.Results.Tol;
    log = p.Results.Log;

    res1 = evalin('caller', str1);
    res2 = evalin('caller', str2);

    nrm = norm(res1 - res2);
    msg = ['norm( [ ', str1, ' ] - [ ', str2, ' ] ) < ', num2str(TOL)];
    if ~isempty(log)
        log_and_print(log, 'Checking that %s...\n', msg);
        log_and_print(log, ['( = ', num2str(nrm), ')\n']);
    end
    if nrm > TOL
        if isempty(log)
            warning(msg)
        else
            log_and_print(log, 'Failed!\r\n');
        end
        pass = false;
    else
        if ~isempty(log)
            log_and_print(log, 'Ok.\r\n');
        end
        pass = true;
    end
    
    if nargout > 1
        nrm_out = nrm;
    end

end

