%-------------------------------------
% BEGIN: function infHorizonDae.m
%-------------------------------------
function dae = infHorizonDae(sol)

tau = sol.time;
y = sol.state;
u = sol.control;

dt_dtau = 2./(1-tau).^2;

ydot = y.*log(y) + y.*u;

dae = [dt_dtau.*ydot];
%-----------------------------------
% END: function infHorizonDae.m
%-----------------------------------

