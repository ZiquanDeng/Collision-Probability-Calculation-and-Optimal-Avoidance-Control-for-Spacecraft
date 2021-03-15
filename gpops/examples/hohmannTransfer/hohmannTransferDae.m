%--------------------------------------%
% BEGIN: function hohmannTransferDae.m %
%--------------------------------------%
function dae = hohmannTransferDae(sol);

global CONSTANTS

t = sol.time;
x = sol.state;
u = sol.control;
p = sol.parameter;
iphase = sol.phase;
r = x(:,1:3);
rad = sqrt(sum(r.*r,2));
er = [r(:,1)./rad r(:,2)./rad r(:,3)./rad];
v = x(:,4:6);
rdot = v;
vdot =-CONSTANTS.mu*[er(:,1)./rad.^2 er(:,2)./rad.^2 er(:,3)./rad.^2];

dae = [rdot vdot];

%-------------------------------------%
% END: function hohmannTransferDae.m  %
%-------------------------------------%
