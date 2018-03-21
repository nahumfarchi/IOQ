function [xHat,wHat,nv,nf,nl] = SEA_det(m,xMax,y,H,cplx)
%
% [xHat,wHat,nv,nf,nl] = SEA_det(m,xMax,y,H,cplx)
%
% Given inputs -
% 
% m    : the problem dimension, i.e., the number of symbols to decode,
% xMax : a parameter specifying the admissible solution space, e.g.,
%        each element of xHat can take values in {-xMax+1,..,-1,0,1,..,xMax},
% y    : a real or complex target signal (column) vector,
% H    : a real or complex linear transform matrix, e.g., a MIMO channel, and
% cplx : flag specifying if xHat should be (0) real- or (1) complex-valued.
% 
% SEA_det returns -
% 
% xHat : a possibly suboptimal solution to   argmin    |y - H*x|^2,
%                                          x \in Z_*^m 
%        where Z_*^m is the set of integers such that -xMax < Z_* <= xMax,
% wHat : the squared distance |y - H*xHat|^2, i.e., the solution node weight,
% nv   : the number of nodes expanded by the detector,
% nf   : the (approximate) number of floating points operations performed
%        (excluding pre-processing, e.g., QR factorization),
% nl   : the number of leaf nodes visited by the detector.
%   
% The real-valued version of this depth-first stack-based sequential decoding
% algorithm uses a decoding tree of m+1 levels where each node has 2*xMax
% children.  At each stage, the node under consideration is expanded if its
% weight is less than the squared distance to nearest currently known lattice
% point.  This distance threshold is initially set to infinity.
%
% Because it is a depth-first traversal, we expand a node by computing its
% first child.  If it is a leaf node, clearly it cannot be further expanded.
% In this case, we will have found a closer lattice point than that previously
% known. Therefore we can adaptively reduce the distance threshold to reflect
% this new discovery.
%
% If the weight of the node under consideration is larger than this distance
% threshold, then the current search path is terminated because it cannot
% possibly lead to a closest lattice point. Upon path termination, the next
% node to be considered is the next sibling of its parent.
%
% If xMax == 0, we do not apply (rectangular, or any) boundary control. In
% other words, SEA_det behaves as a lattice decoder.  
%
% For more sophisticated operation, xMax may also be a vector of length m.
% Then each node at the beginning of stage j in the tree, where the root node
% is at the beginning of stage m and the leaf nodes are found at the end of
% stage 1, has 2*xMax(j) children.  Equivalently, symbol xHat(j) is drawn from
% {-xMax(j)+1,..,-1,0,1,..,xMax(j)}.
%
% If cplx == 1, we consider a tree of 2*m+1 levels with each node still having
% 2*xMax children.  In addition, either xMax should be a complex-valued vector,
% or imag(xMax) will be taken to be equal to real(xMax), i.e., a square QAM
% constellation will be assumed by default.
% 
% If either lattice reduction assistance or MMSE pre-processing are desired,
% these operations should also be applied in advance of calling SEA_det.
%
% Notes:
%
% - Size(H,1) is expected to be equal to length(y).
% - Size(H,1) is expected to be greater than or equal to size(H,2).
% - Size(H,2) is expected to be equal to m.
% - Length(xMax) is expected to be equal to either 1 or m.
% - The solution xHat is an integer vector; be sure to apply any necessary
%   scaling prior to calling SEA_det so that this solution is appropriate.
%
% Example:
%
%   m    = 4;                        % 4x4 MIMO system
%   xMax = 2;                        % 16-QAM modulation
%   cplx = 1;
%   y    = randn(m,1)+i(randn(m,1);  % generate random complex target vector
%   H    = randn(m,m)+i*randn(m,m);  % generate random complex channel 
%   [xHat,wHat,nv,nf,nl] = SEA_det(m,xMax,y,H,cplx);
%
% Version 1.2, copyright 2006 by Karen Su.
%
% Please send comments and bug reports to karen.su@utoronto.ca.
%

% Version 1.1
%   Bug fix: now correctly handles xMax == 0 case.
% Version 1.2
%   Bug fix: Now correctly counts number of nodes visited; v1.1 did not count
%            first expansion and so under-reported nv by exactly 1.

global w;
global z;
global cn;
global zp;
global nf;

xHat = [];
wHat = Inf;
nv   = 0;
nf   = 0;
nl   = 0;

%
% A small number
%
EPS = 0.001;

%
% Basic error checking
%
[n,mChk] = size(H);
if n ~= length(y) || mChk ~= m
  fprintf('Error: Decoding failed, invalid dimensions for H.\n');
  return;
end
if n < m
  fprintf('Error: Cannot solve under-determined problem, n (%i) < m (%i).\n', n, m);
  return;
end
if length(xMax) == 1
  xMax = xMax*ones(m,1);
elseif size(xMax,2) == m
  xMax = xMax.';
elseif size(xMax,1) ~= m
  fprintf('Error: Decoding failed, invalid dimensions for xMax.\n');
  return;
end

%
% Complex -> real conversion if necessary
%
if cplx == 1
  m = 2*m;
  n = 2*n;
  y = [real(y);imag(y)];
  H = [[real(H),-imag(H)];[imag(H),real(H)]];
  if isreal(xMax)
    xMax = [xMax;xMax];
  else
    xMax = [real(xMax);imag(xMax)];
  end
end

%
% Special handling to allow xMax == 0
%
noBCind       = find(xMax==0);
xMax(noBCind) = Inf;

%
% Reduce problem to square if H is overdetermined; also factorize code
% generator
%
[Hq,Hr] = qr(H);
wRoot   = 0;
yR      = Hq'*y;
if m < n
  wRoot = norm(yR(m+1:n))^2;
  yR    = yR(1:m);
end

%
% A basic node data structure consists of
%   w   : node weight
%   pyR : parent's residual target vector, elements 1:dim  == py(1:dim,dim)
%   zA  : accumulated symbol decision vector, length m-dim == z(dim+1:dim)
%   dim : node/residual dimension, root (m+1) to leaf (1)
%   cn  : sibling ordinal
%   pw  : parent node weight                               == w(dim+1)
%   ut  : (scalar) target for siblings at this level       == py(dim+1,dim)
%   zp  : symbol decision of previous sibling             
%
% In addition, the following components enable us to effect more sophisticated,
% albeit still justified rectangular, boundary control
%   lb  : lower bound of admissible range for siblings at this level == lb(dim)
%   ub  : upper bound of admissible range for siblings at this level == ub(dim)
%
% The algorithm requires only a fixed-size block of memory; we initialize it
% with the root node.  
%
w   = zeros(1,m+1);
py  = zeros(m,m);
z   = zeros(1,m);
dim = m;
cn  = zeros(1,m);
zp  = zeros(1,m);

lb  = -xMax.'+1;
ub  = xMax.';

dimp    = dim+1;
w(dimp) = wRoot;

%
% Compute first child of root node 
%
py(1:dim,dim) = yR;                              % Store parent (root) residual target.
compFirstChild(dim,dimp,lb(dim),ub(dim),Hr(dim,dim),py(dim,dim));
nv = nv + 1;

backtrack = 0;
while dim <= m
  if w(dim) < wHat                     % If the node under consideration has a smaller weight
    if dim == 1                        % and it is a leaf node, then update xHat, wHat 
      wHat = w(dim);
      xHat = z;
      nl   = nl + 1;

      backtrack = 1;                   % To backtrack, we have to compute the next sibling
      dim  = dimp;                     % of the current node's parent
      dimp = dimp+1;
      if cn(dim) <= ub(dim)-lb(dim)    % Check that there is at least one more sibling                 
        backtrack = 0;
        compNextSibling(dim,dimp,lb(dim),ub(dim),Hr(dim,dim),py(dim,dim));
      end
      nf = nf + 1;

    %
    % If the node under consideration has a sufficiently small weight and is
    % not a leaf node:
    % - then either we are backtracking (backtrack == 1, and so compute the
    %   next sibling of its parent)
    %
    elseif backtrack && cn(dim) > ub(dim)-lb(dim) % No more siblings
      %backtrack = 1;
      dim = dimp;
      if dim <= m                      % Check that we have not backtracked to root
        dimp = dimp+1; 
        if cn(dim) <= ub(dim)-lb(dim)  % Check that there is at least one more sibling
          backtrack = 0;
          compNextSibling(dim,dimp,lb(dim),ub(dim),Hr(dim,dim),py(dim,dim));
        end
        nf = nf + 1;
      end
  
    %
    % - or we are trying an alternate path (backtrack == 0, and so compute
    %   the first child)
    %
    else
      %backtrack = 0;
      nv   = nv + 1;
      dimp = dim;
      dim  = dim-1;
      py(1:dim,dim) = py(1:dim,dimp)-Hr(1:dim,dimp)*z(dimp);  % Compute current node residual target.
      nf            = nf + 2*dim;
      compFirstChild(dim,dimp,lb(dim),ub(dim),Hr(dim,dim),py(dim,dim));
    end

  else                                 % If the node under consideration has too large a weight
    dim       = dimp;                  % then we backtrack (if possible).  To do so, we have to
    backtrack = 1;                     % compute the next sibling of the current node's parent.

    if dim <= m                        % Check that we have not backtracked to root
      dimp = dimp+1;
      if cn(dim) <= ub(dim)-lb(dim)    % Check that there is at least one more sibling
        backtrack = 0;
        compNextSibling(dim,dimp,lb(dim),ub(dim),Hr(dim,dim),py(dim,dim));
      end
      nf = nf + 1;
    end
  end
  nf = nf + 1;
end

xHat = xHat.';

%   
% Undo any complex -> real transformation that may have been done before
%   
if cplx == 1
  m = m/2;
  xHat = xHat(1:m)+i*xHat(m+1:2*m);
end

function compFirstChild(dim,dimp,lbd,ubd,Hrdd,pydd)
global w;
global z;
global cn;
global zp;
global nf;

z0 = pydd/Hrdd;                              % The first child is found by
z1 = round(z0);                              % quantization and symbol-dependent
z1 = min([z1 ubd]);                          % boundary control.
z1 = max([z1 lbd]);

z(dim)  = z1;                                % Store first child decision
zp(dim) = z0;                                % and pre-quantization data;
w(dim)  = w(dimp)+(pydd-Hrdd*z(dim))^2;      % compute its weight/cost.
cn(dim) = 1;                                 % First child ordinal has 1.
nf      = nf + 6;

function compNextSibling(dim,dimp,lbd,ubd,Hrdd,pydd)
global w;
global z;
global cn;
global zp;
global nf;

zn   = z(dim)-cn(dim)*sign(z(dim)-zp(dim));  % The next sibling symbol decision can be
if zn > ubd                                  % derived from the current and previous ones,
  zn = ubd-cn(dim);                          % followed by boundary control. This procedure
elseif zn < lbd                              % is dependent on the child ordinal and also the
  zn = lbd+cn(dim);                          % symbol-dependent lower and upper bounds.
end
zp(dim) = z(dim);                            % Save previous sibling's symbol decision and
z(dim)  = zn;                                % store next sibling decision.
cn(dim) = cn(dim)+1;                         % Siblings have one larger child ordinal.

w(dim) = w(dimp)+(pydd-Hrdd*zn)^2;           % Compute sibling node weight.
nf     = nf + 10;