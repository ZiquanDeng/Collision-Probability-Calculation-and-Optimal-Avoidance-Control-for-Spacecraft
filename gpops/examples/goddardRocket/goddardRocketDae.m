function dae = goddardRocketDae(sol)

global CONSTANTS

t = sol.time;
x = sol.state;
u = sol.control;
p = sol.parameter;
h = x(:,1);
v = x(:,2);
m = x(:,3);
T = sol.control;
D = CONSTANTS.dragk.*(v.^2).*exp(-h/CONSTANTS.H);
hdot = v;
vdot = (T-D)./m-CONSTANTS.g0*ones(size(t));
mdot = -T./CONSTANTS.c;
if sol.phase==2,
    voverc = v/CONSTANTS.c;
    xmg = m*CONSTANTS.g0;
    term1 = (CONSTANTS.c^2).*(ones(size(t))+voverc)./(CONSTANTS.H*CONSTANTS.g0)-ones(size(t))-2./voverc;
    term2 = xmg./(ones(size(t))+4./voverc+2./(voverc.^2));
    path = T-D-xmg-term1.*term2;
    dae = [hdot vdot mdot path];
else
    dae = [hdot vdot mdot];
end;

%-------------------------------%
% End File:  goddardRocketDae.m %
%-------------------------------%