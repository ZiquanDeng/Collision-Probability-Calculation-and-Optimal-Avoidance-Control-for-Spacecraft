%--------------------------------------
% BEGIN: function infHorizonCost.m
%--------------------------------------
function [Mayer,Lagrange]=infHorizonCost(sol)

tau  = sol.time;
y = sol.state;
u = sol.control;

dt_dtau = 2./(1-tau).^2;

Mayer = 0;
Lagrange = 0.5 * dt_dtau .* ( log(y).^2 + u.^2 );

%--------------------------------------
% END: function infHorizonCost.m
%--------------------------------------