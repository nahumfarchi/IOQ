
function z_hat = babai(R, y)
%%
%   compute the Babai estimation
%   find a sub-optimal solution for min_z ||R*z-y||_2
%   B - lattice basis
%   y - a real vector of n-by-1
%   z_hat - resulting integer vector
%

    %R = LLL_reduction(R);
    %R = chol(B);
    n=length(y);
    z_hat=zeros(n,1);
    z_hat(n)=round(y(n)./R(n,n));
    for k=n-1:-1:1
        par=R(k,k+1:n)*z_hat(k+1:n);
        ck=(y(k)-par)./R(k,k);
        z_hat(k)=round(ck);
    end

%     w = w(:);
%     n = size(V, 2);
%     m = size(V, 1);
%     B = LLL_reduction(V);
%     
%     [B_tilde, ~] = qr(B);
%     R = zeros(m, n);
%     for i = 1:m
%         bi = B(:, i);
%         for j = 1:n
%             if i <= j
%                 bi_tilde = B_tilde(:, i);
%                 bj = B(:, j);
%                 muij = dot(bi_tilde, bj) / norm(bi)^2;
%                 R(i,j) = norm(bi) / muij;
%             end
%         end
%     end
% 
%     y = w;
%     n=length(y);
%     z_hat=zeros(n,1);
%     z_hat(n)=round(y(n)./R(n,n));
%     for k=n-1:-1:1
%         par=R(k,k+1:n)*z_hat(k+1:n);
%         ck=(y(k)-par)./R(k,k);
%         z_hat(k)=round(ck);
%     end
%     x = z_hat;

%     
%     t = R \ w;
%     s = -t;
%     for i = 1:n
%         b_tilde = B_tilde(:, n-i+1);
%         c = round(dot(s, b_tilde) / norm(b_tilde)^2);
%         s = s - c*B(:, n-i+1);
%     end
%     x = R*s;

    
    
    %w = B_tilde' * w;
%     b = w;
%     for j = n:-1:1
%        bj_tilde = B_tilde(:, j);
%        bj = B(:, j);
%        cj = round(dot(b, bj_tilde) / dot(bj_tilde, bj_tilde));
%        b = b - cj*bj;
%     end
%     x = w - b;
    %x = B_tilde * x;
    
    %R = LLL_reduction(B);
    %[~, R] = qr(B_reduced);
    %y = R \ y;
    
    %n=length(y);
    %z_hat=zeros(n,1);
    %z_hat(n)=round(y(n)./R(n,n));
    %for k=n-1:-1:1
    %    par=R(k,k+1:n)*z_hat(k+1:n);
    %    ck=(y(k)-par)./R(k,k);
    %    z_hat(k)=round(ck);
    %end

end