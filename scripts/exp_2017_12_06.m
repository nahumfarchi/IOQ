%Program for finding the inverse of a given matrix in two steps
%Opening of the matrix to be inverted
%A = xlsread('MatrizA.xls');
%Alternative creation of a random matrix
% prompt = {'Enter "Y" for a random matrix:'};
% dlg_title = 'Input for random matrix';
% num_lines = 1;
% defaultanswer = {'N'};
% answer = inputdlg(prompt,dlg_title,num_lines,defaultanswer);
% resp= char(answer); 
% switch (resp)
%    case 'Y'
%       SM= uint16(input('Desired size of the random matrix: '));
%       A=rand(SM);
%       if SM <= 26 
%         seq= 1;
%         NL1= SM;
%         NL2= 0;
%         NL3= 0;
%       elseif SM <= 676
%         seq= 2;
%         NL1= idivide(SM,26,'fix');
%         NL2= SM- NL1*26;
%         NL3= 0;
%       elseif SM <= 2000 
%         seq= 3;
%         NL1= idivide(SM,676,'fix');
%         NL2= idivide(SM- NL1*676,26,'fix');
%         NL3= SM- NL1*676- NL2*26;
%       else
%         error('The size of the matrix is over 256');
%       end
%       ch1= fSmToExcel(NL1);
%       ch2= fSmToExcel(NL2);
%       ch3= fSmToExcel(NL3);
%       str = int2str(SM);
%       s = ['A1:',ch1,ch2,ch3,str];
%       xlswrite('MatrizA.xls',A,s);
%    case 'N'
%       disp('The imported matrix is used')
%    otherwise
%       disp('The imported matrix is used')
% end
%A = magic(4);

% m = Mesh('../data/bucky180.off');
% [d0, ~] = get_exterior_derivatives(m);
% A = d0' * d0;
% 
% mo = size(A,1);
% disp('Matrix size: ');
% disp(mo);
% B=eye(mo);
% C= [A B]; % creating the amplied matrix
% % Now the rough inverse is calculated by Gauss-Jordan method
% %errperm=input('Maximal permisible absolute error: ');
% errperm = 100;
% tic
% for i= 1:(mo-1);
%     for k=i+1:mo;
%         fac=C(k,i)/C(i,i);
%         for l= 1:2*mo;
%             C(k,l)=C(k,l)-fac*C(i,l);
%         end
%     end
% end
% X = C;
% l=2*mo;
% for l=1:mo;
%    X(:, 1) = [];
% end
% for j=1:mo;
%     i=mo;
%     while i~=0
%         for k=i+1:mo;
%             if (i~=mo)
%                 X(i,j)=X(i,j)-C(i,k)*X(k,j);
%             end;
%         end
%         X(i,j)=X(i,j)/C(i,i);
%         i=i-1;
%     end
% end
% %The matrix X is the rough inverse of A; now starts the refination
% RMAX= 1;
% it= 0;
% while RMAX> errperm
%     B=A*X;
%     R= 2*eye(mo)-B;
%     RMAX= abs(1-norm(R,1));
%     XC=X*R;
%     X=XC;
%     it= it+1;
% end
% TC= toc;
% Nrd=norm(A,1);
% Nri=norm(X,1);
% cond= Nrd*Nri;
% disp('Inverse matrix:  ');
% disp(X);
% disp('Number of iterations:  ');
% disp(it);
% disp('Error of the inverse matrix:  ');
% disp(RMAX);
% disp('Time consumption for running (s):  ');
% disp(TC);
% disp('Conditioning of matrix A:  ');
% disp(cond);

gd = gpuDevice();
m = Mesh('../data/sphere_s4.off');
nv = m.nV
[d0, ~] = get_exterior_derivatives(m);
L = d0' * d0;
Lfull = full(L+1/nv);
disp('starting...')
%[R, p] = chol(L);
tic; 
[R,p]=chol(L+1/nv); 
Rinv = inv(R);
Rinv = Rinv' * Rinv;
%Rinv = inv(R)
toc;
fprintf('err = %g\n', norm(full(speye(nv) - Rinv*Lfull)));

tic
Linv = sparseinv(L);
toc
fprintf('err = %g\n', norm(full(speye(nv) - Linv*L)));
%tic
%[R,p]=chol(L);
%Rinv = inv(R);
%Rinv = Rinv' * Rinv;
%toc
%assert(p == 0)
%timeit(@() [R,p]=chol(L))

% chol mex inv
tic
invChol_mex(full(L+1/nv));
toc

% GPU

% chol
tic; 
[Rgpu, pgpu] = chol(gpuArray(full(L+1/nv))); 
Rgpu = inv(Rgpu);
Rgpu = Rinvgpu' * Rinvgpu;
wait(gd); toc
clear Rgpu;
clear pgpu;

% backslash s4: -
tic
Lgpuinv = gpuArray(full(L+1/nv)) \ eye(nv);
wait(gd); toc
clear Lgpuinv

% inv s4: -
tic
Lgpuinv = inv(gpuArray(full(L+1/nv)));
wait(gd); toc
clear Lgpuinv

% block inv s4: 177
tic;
block_inv_gpu(full(L+1/nv), 8000, true);
wait(gd); toc
%assert(p == 0)
%gputimeit(@() chol(Lgpu))