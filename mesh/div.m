function [divop] = div(m, gradop)
   if nargin < 2
       [gradop, ~] = grad(m);
   end
   divop = -m.Gv_inv*gradop'*m.Gf;
end

