function [C,G] = gpopsuserfunCS(x)
%------------------------------------------------------------------%
% User function when complex step differentiation is used          %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

global mysetup

%--------------------------------%
% Unscale the decision variables %
%--------------------------------%
% y = x./mysetup.column_scales;
y = (x-mysetup.column_shifts)./mysetup.column_scales;
[C,J] = gpopsObjandConsComplex(y);
%-----------------------------------------------------------------%
% Compute functions based on unscaled value of decision variables %
%-----------------------------------------------------------------%
C(2:mysetup.numnonlin+1) = mysetup.DF*C(2:mysetup.numnonlin+1);
J = J*mysetup.invDx;
%----------------------%
% Unscale the Jacobian %
%----------------------%
J(2:mysetup.numnonlin+1,:) = mysetup.DF*J(2:mysetup.numnonlin+1,:);
G = snfindG(mysetup.iGfun,mysetup.jGvar,J);

