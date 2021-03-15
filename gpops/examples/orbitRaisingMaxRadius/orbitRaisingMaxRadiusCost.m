function [Mayer,Lagrange]=orbitRaisingMaxRadiusCost(solcost);

t0 = solcost.initial.time;
x0 = solcost.initial.state;
tf = solcost.terminal.time;
xf = solcost.terminal.state;
t  = solcost.time;
x  = solcost.state;
u  = solcost.control;
p  = solcost.parameter;

Mayer = -xf(1);
Lagrange = zeros(size(t)); 
