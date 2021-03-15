%-------------------------------------------%
% Begin Function:  brysonMinimumClimbCost.m %
%-------------------------------------------%
function [Mayer,Lagrange] = brysonMinimumClimbCost(solcost);

tf = solcost.terminal.time;
t  = solcost.time;

Mayer = tf;
Lagrange = zeros(size(t));
%-----------------------------------------%
% End Function:  brysonMinimumClimbCost.m %
%-----------------------------------------%
