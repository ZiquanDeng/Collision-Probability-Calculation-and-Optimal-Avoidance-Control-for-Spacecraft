%-------------------------------------%
% BEGIN: function hyperSensitiveDae.m %
%-------------------------------------%
function [dae] = hyperSensitiveDae(sol);

t = sol.time;
x = sol.state;
u = sol.control;
p = sol.parameter;

dae = -x.^3+u;

%-------------------------------------%
% END: function hyperSensitiveDae.m   %
%-------------------------------------%
