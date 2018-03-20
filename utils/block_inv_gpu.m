function Xinv = block_inv_gpu(X, lb, use_block_mult)

    if nargin < 3
        use_block_mult = true;
    end

    % 60s on 57k
    n = size(X, 1);
    n2 = round(n / 2);
    %fprintf("n, n2 = %d, %d\n", n, n2);
    
    %A = X(1:n2,1:n2);
    %B = X(1:n2,n2+1:end);
    %C = X(n2+1:end,1:n2);
    %D = X(n2+1:end,n2+1:end);
    
    if n2 > lb
        Ainv = gpuArray(block_inv_gpu(X(1:n2,1:n2), lb));
        
        if use_block_mult
            E = block_mult_gpu(Ainv, X(1:n2,n2+1:end), lb, 'single');
            Ainv = gather(Ainv);
            E = block_mult_gpu(X(n2+1:end,1:n2), E, lb, 'single');
            E = gpuArray(X(n2+1:end,n2+1:end)) - gpuArray(E);
            E = gather(E);
        else
            E = gather(X(n2+1:end,n2+1:end) - ...
                X(n2+1:end,1:n2)*Ainv*X(1:n2,n2+1:end)); % !
            Ainv = gather(Ainv);
        end
        %E = X(n2+1:end,1:n2)*Ainv*X(1:n2,n2+1:end);
        %E = X(n2+1:end,1:n2)*Ainv*X(1:n2,n2+1:end);
        %E = gather(E);
        
        
        
        %E = gather(X(n2+1:end,n2+1:end) - ...
        %    X(n2+1:end,1:n2)*(gpuArray(X(1:n2,1:n2))\X(1:n2,n2+1:end)));
        %Ainv = gather(Ainv);
        Einv = block_inv_gpu(E, lb);
        
        %Einv = gpuArray(Einv);
        Ainv = gpuArray(Ainv);
        %X = gpuArray(X);
    else
        %X = gpuArray(X);
        Ainv = inv(gpuArray(X(1:n2,1:n2)));
        Einv = gather(inv(X(n2+1:end,n2+1:end) - X(n2+1:end,1:n2) * ...
            Ainv * X(1:n2,n2+1:end)));
        %Einv = gather(inv(X(n2+1:end,n2+1:end) - X(n2+1:end,1:n2) * ...
        %    (gpuArray(X(1:n2,1:n2)) \ X(1:n2,n2+1:end))));
        %Ainv = gather(Ainv);
    end
    
    % Here only Ainv is on the gpu
    
    Xinv = zeros(n, n, 'single');
    if n2 < lb || ~use_block_mult
        %Ainv = gpuArray(Ainv);
        %Einv = gpuArray(Einv);
        Xinv(1:n2,1:n2) = gather(Ainv + ...
            Ainv * X(1:n2,n2+1:end) * Einv * X(n2+1:end,1:n2) * Ainv);
        Xinv(1:n2,n2+1:end) = gather(-Ainv * X(1:n2,n2+1:end) * Einv);
        Xinv(n2+1:end,1:n2) = gather(-Einv * X(n2+1:end,1:n2) * Ainv);
    else
        
        tmp = block_mult_gpu(Ainv, X(1:n2,n2+1:end), lb, 'single');
        tmp = block_mult_gpu(tmp, Einv, lb, 'single');
        tmp = block_mult_gpu(tmp, X(n2+1:end,1:n2), lb, 'single');
        tmp = block_mult_gpu(tmp, Ainv, lb, 'single');
        Xinv(1:n2,1:n2) = gather(Ainv + gpuArray(tmp));
        
        tmp = block_mult_gpu(-Ainv, X(1:n2,n2+1:end), lb, 'single');
        Xinv(1:n2,n2+1:end) = block_mult_gpu(tmp, Einv, lb, 'single');
        
        tmp = block_mult_gpu(-Einv, X(n2+1:end,1:n2), lb, 'single');
        Xinv(n2+1:end,1:n2) = block_mult_gpu(tmp, Ainv, lb, 'single');
    end
    clear Ainv;
    Xinv(n2+1:end,n2+1:end) = gather(Einv);
    
    %fprintf("exiting %d, %d\n", n, n2);


% backslash
%     n = size(X, 1);
%     n2 = round(n / 2);
%     
%     %A = X(1:n2,1:n2);
%     %B = X(1:n2,n2+1:end);
%     %C = X(n2+1:end,1:n2);
%     %D = X(n2+1:end,n2+1:end);
%     
%     if n2 > lb
%         Ainv = (block_inv_gpu(X(1:n2,1:n2), lb));
%         %E = gather(X(n2+1:end,n2+1:end) - X(n2+1:end,1:n2)*Ainv*X(1:n2,n2+1:end));
%         E = gather(X(n2+1:end,n2+1:end) - ...
%             X(n2+1:end,1:n2)*(gpuArray(X(1:n2,1:n2))\X(1:n2,n2+1:end)));
%         %Ainv = gather(Ainv);
%         Einv = block_inv_gpu(E, lb);
%         
%         %Einv = gpuArray(Einv);
%         %Ainv = gpuArray(Ainv);
%         %X = gpuArray(X);
%     else
%         %X = gpuArray(X);
%         Ainv = inv(gpuArray(X(1:n2,1:n2)));
%         Einv = gather(inv(X(n2+1:end,n2+1:end) - X(n2+1:end,1:n2) * ...
%             (X(1:n2,1:n2) \ X(1:n2,n2+1:end))));
%         %Einv = gather(inv(X(n2+1:end,n2+1:end) - X(n2+1:end,1:n2) * ...
%         %    (gpuArray(X(1:n2,1:n2)) \ X(1:n2,n2+1:end))));
%         %Ainv = gather(Ainv);
%     end
%     
%     % Here only Ainv is on the gpu
%     
%     Xinv = zeros(n, n, 'single');
%     %Ainv = gpuArray(Ainv);
%     %Einv = gpuArray(Einv);
%     Xinv(1:n2,1:n2) = gather(X(1:n2,1:n2) \ (eye(n2) + ...
%         Ainv * X(1:n2,n2+1:end) * Einv * X(n2+1:end,1:n2) * Ainv));
%     Xinv(1:n2,n2+1:end) = gather((-X(1:n2,1:n2) \ X(1:n2,n2+1:end)) * Einv);
%     Xinv(n2+1:end,1:n2) = gather(-Einv * (X(n2+1:end,1:n2) / X(1:n2,1:n2)));
%     clear Ainv;
%     Xinv(n2+1:end,n2+1:end) = gather(Einv);


% assert(size(X,1) == size(X,2))
%     n = size(X, 1);
%     n2 = round(n/2);
%     
%     %A = X(1:n2,1:n2);
%     %B = X(1:n2,n2+1:end);
%     %C = X(n2+1:end,1:n2);
%     %D = X(n2+1:end,n2+1:end);
%     if n2 > lb
%         %Xinv = zeros(n, n, 'gpuArray');
%         %Xinv(1:n2,1:n2) = block_inv(X(1:n2,1:n2), lb);
%         %Xinv(1:n2,n2+1:end) = block_inv(X(1:n2,n2+1:end), lb);
%         %Xinv(n2+1:end,1:n2) = block_inv(X(n2+1:end,1:n2), lb);
%         %Xinv(n2+1:end,n2+1:end) = block_inv(X(n2+1:end,n2+1:end), lb);
%         Ainv = block_inv_gpu(X(1:n2,1:n2), lb);
%         
%         %B = X(1:n2,n2+1:end);
%         %C = X(n2+1:end,1:n2);
%         %D = X(n2+1:end,n2+1:end);
%         %Einv = block_inv(D-C*Ainv*B, lb);
%         X = gpuArray(X);
%         Einv = block_inv_gpu(X(n2+1:end,n2+1:end)-X(n2+1:end,1:n2)*Ainv*X(1:n2,n2+1:end), lb);
%     else
%         X = gpuArray(X);
%         Ainv = inv(X(1:n2,1:n2));
%         
%         %Einv = inv(D-C*Ainv*B);    
%         Einv = inv(X(n2+1:end,n2+1:end) - X(n2+1:end,1:n2)*Ainv*X(1:n2,n2+1:end));
%     end
%     
%     tmp = Ainv*X(1:n2,n2+1:end)*Einv;
%     %clear A;
%     %clear B;
%     %clear C;
%     %clear D;
%     %X = gather(X);
%     %Xinv = zeros(n, n, 'gpuArray');
%     
% %     % suspect
% %     X(1:n2, 1:n2) = Ainv + tmp*X(n2+1:end,1:n2)*Ainv;
% %     X(n2+1:end, 1:n2) = -Einv*X(n2+1:end,1:n2)*Ainv;
% %     clear Ainv;
% %     X(1:n2, n2+1:end) = -tmp;
% %     clear tmp;
% %     X(n2+1:end, n2+1:end) = Einv;
%     
%     %X11 = Ainv + tmp*X(n2+1:end,1:n2)*Ainv;
%     %X21 = -Einv*X(n2+1:end,1:n2)*Ainv;
%     %clear Ainv;
%     %X = gather(X);
%     %X12 = -tmp;
%     %X22 = Einv;
%     X = [Ainv + tmp*X(n2+1:end,1:n2)*Ainv , -tmp;
%                -Einv*X(n2+1:end,1:n2)*Ainv  , Einv];

%     % 48
%     n = size(X, 1);
%     n2 = round(n / 2);
%     
%     %A = X(1:n2,1:n2);
%     %B = X(1:n2,n2+1:end);
%     %C = X(n2+1:end,1:n2);
%     %D = X(n2+1:end,n2+1:end);
%     
%     if n2 > lb
%         Ainv = gpuArray(block_inv_gpu(X(1:n2,1:n2), lb));
%         E = gather(X(n2+1:end,n2+1:end) - X(n2+1:end,1:n2)*Ainv*X(1:n2,n2+1:end));
%         %E = gather(X(n2+1:end,n2+1:end) - ...
%         %    X(n2+1:end,1:n2)*(gpuArray(X(1:n2,1:n2))\X(1:n2,n2+1:end)));
%         Ainv = gather(Ainv);
%         Einv = block_inv_gpu(E, lb);
%         
%         %Einv = gpuArray(Einv);
%         Ainv = gpuArray(Ainv);
%         %X = gpuArray(X);
%     else
%         %X = gpuArray(X);
%         Ainv = inv(gpuArray(X(1:n2,1:n2)));
%         Einv = (inv(X(n2+1:end,n2+1:end) - X(n2+1:end,1:n2) * ...
%             Ainv * X(1:n2,n2+1:end)));
%         %Einv = gather(inv(X(n2+1:end,n2+1:end) - X(n2+1:end,1:n2) * ...
%         %    (gpuArray(X(1:n2,1:n2)) \ X(1:n2,n2+1:end))));
%         %Ainv = gather(Ainv);
%     end
%     
%     % Here both Ainv and Einv are on the gpu
%     
%     Xinv = zeros(n, n, 'single');
%     Xinv(n2+1:end,n2+1:end) = gather(Einv); % D
%     Xinv(1:n2,n2+1:end) = gather(-Ainv*X(1:n2,n2+1:end)*Einv); % B
%     Einv = Einv*X(n2+1:end,1:n2)*Ainv;
%     Xinv(n2+1:end,1:n2) = gather(-Einv); % C
%     Xinv(1:n2,1:n2) = gather(Ainv + Ainv*X(1:n2,n2+1:end)*Einv);
    
    %Ainv = gpuArray(Ainv);
    %Einv = gpuArray(Einv);
%     Xinv(1:n2,1:n2) = gather(Ainv + Ainv * X(1:n2,n2+1:end) * Einv * X(n2+1:end,1:n2) * Ainv);
%     Xinv(1:n2,n2+1:end) = gather(-Ainv * X(1:n2,n2+1:end) * Einv);
%     Xinv(n2+1:end,1:n2) = gather(-Einv * X(n2+1:end,1:n2) * Ainv);
%     clear Ainv;
%     Xinv(n2+1:end,n2+1:end) = gather(Einv);