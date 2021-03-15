function daeout = moonlanderDae(soldae);

global CONSTANTS

t = soldae.time;
x = soldae.state;
w = soldae.control;
p = soldae.parameter;

h = x(:,1);
v = x(:,2);
u = soldae.control;
hdot = v;
vdot = -CONSTANTS.g+u;
daeout = [hdot vdot];

%-------------------------------%
% END: function moonlanderDae.m %
%-------------------------------%
