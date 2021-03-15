%---------------------------------------%
% BEGIN: function brachistochroneCost.m %
%---------------------------------------%
function [Mayer,Lagrange]=brachistochroneCost(sol);

tf = sol.terminal.time;
t  = sol.time;

Mayer = tf;
Lagrange = zeros(size(t));

%-------------------------------------%
% END: function brachistochroneCost.m %
%-------------------------------------%
