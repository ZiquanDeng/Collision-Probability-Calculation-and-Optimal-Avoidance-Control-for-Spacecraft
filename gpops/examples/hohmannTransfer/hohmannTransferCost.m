%---------------------------------------%
% BEGIN: function hohmannTransferCost.m %
%---------------------------------------%
function [Mayer,Lagrange] = hohmannTransferCost(sol);

global CONSTANTS
t = sol.time;
t0 = sol.initial.time;
x0 = sol.initial.state;
tf = sol.terminal.time;
xf = sol.terminal.state;
p = sol.parameter;
Mayer = p(1:3).'*p(1:3)+p(4:6).'*p(4:6);
Lagrange = zeros(size(t));

%-------------------------------------%
% END: function hohmannTransferCost.m %
%-------------------------------------%
