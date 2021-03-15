function [Mayer,Lagrange]=linearTangentCost(sol);

t0 = sol.initial.time;
x0 = sol.initial.state;
tf = sol.terminal.time;
xf = sol.terminal.state;
t  = sol.time;
x  = sol.state;
u  = sol.control;
p  = sol.parameter;

Mayer = tf;
Lagrange = zeros(size(t));

%-----------------------------------%
% END: function linearTangentCost.m %
%-----------------------------------%
