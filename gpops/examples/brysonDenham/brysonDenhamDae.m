%-----------------------------------%
% BEGIN: function brysonDenhamDae.m %
%-----------------------------------%
function [dae] = brysonDenhamDae(sol);

t = sol.time;
x = sol.state;
u = sol.control;
x1dot = x(:,2);
x2dot = u;
x3dot = u.^2/2;
dae = [x1dot x2dot x3dot];

%-----------------------------------%
% END: function brysonDenhamDae.m   %
%-----------------------------------%
