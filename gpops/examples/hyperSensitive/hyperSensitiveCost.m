%--------------------------------------%
% BEGIN: function hyperSensitiveCost.m %
%--------------------------------------%
function [Mayer,Lagrange] = hypersensitiveCost(sol);

t0 = sol.initial.time;
x0 = sol.initial.state;
tf = sol.terminal.time;
xf = sol.terminal.state;
t  = sol.time;
x  = sol.state;
u  = sol.control;
p  = sol.parameter;
Mayer = zeros(size(t0));
Lagrange = 0.5*(x.^2+u.^2);

%--------------------------------------%
% END: function hyperSensitiveCost.m   %
%--------------------------------------%
