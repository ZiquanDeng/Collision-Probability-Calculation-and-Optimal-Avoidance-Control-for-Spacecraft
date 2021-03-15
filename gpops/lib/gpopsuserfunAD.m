function [C,G] = gpopsuserfunAD(x)
%------------------------------------------------------------------%
% User function when automatic differentiation is used             %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

global mysetup

%--------------------------------%
% Unscale the decision variables %
%--------------------------------%
y = ad(x);
y = (x-mysetup.column_shifts)./mysetup.column_scales;
%-----------------------------------------------------------------%
% Compute functions based on unscaled value  of decision variables %
%-----------------------------------------------------------------%
CC = gpopsObjandCons(y);
%----------------------------------------------------------------%
% Get the values of the constraint vector and objective function %
%----------------------------------------------------------------%
C = getvalue(CC);
%-----------------------%
% Scale the constraints %
%-----------------------%
C(2:mysetup.numnonlin+1) = mysetup.DF*C(2:mysetup.numnonlin+1);
%------------------------------------------------------------%
% Get the Jacobian of the constraints and objective function %
%------------------------------------------------------------%
J = getderivative(CC)*mysetup.invDx;
%----------------------%
% Unscale the Jacobian %
%----------------------%
J(2:mysetup.numnonlin+1,:) = mysetup.DF*J(2:mysetup.numnonlin+1,:);
G = snfindG(mysetup.iGfun,mysetup.jGvar,J);
