function [Mayer,Lagrange] = tumorAntiAngiogenesisCost(solcost)

t0 = solcost.initial.time;
x0 = solcost.initial.state;
tf = solcost.terminal.time;
xf = solcost.terminal.state;
t  = solcost.time;
x  = solcost.state;
u  = solcost.control;
p  = solcost.parameter;

pf = xf(1);

Mayer = pf;
Lagrange = zeros(size(t));
