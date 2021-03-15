function [Mayer,Lagrange]=moonlanderCost(solcost);

t0 = solcost.initial.time;
x0 = solcost.initial.state;
tf = solcost.terminal.time;
xf = solcost.terminal.state;
t  = solcost.time;
x  = solcost.state;
w  = solcost.control;
p  = solcost.parameter;

u = solcost.control;
Mayer = zeros(size(t0));
Lagrange = u;

%--------------------------------%
% END: function moonlanderCost.m %
%--------------------------------%
