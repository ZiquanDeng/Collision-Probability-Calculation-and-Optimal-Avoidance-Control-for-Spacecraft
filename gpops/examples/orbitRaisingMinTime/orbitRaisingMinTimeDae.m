function dae = orbitRaisingMinTimeDae(soldae,iphase);

global CONSTANTS

t = soldae.time;
x = soldae.state;
u = soldae.control;
p = soldae.parameter;

r      = x(:,1);
vr     = x(:,2);
vtheta = x(:,3);
m      = x(:,4);
u1     = u(:,1);
u2     = u(:,2);

rdot = vr;
vrdot = (vtheta.^2)./r-CONSTANTS.mu./(r.^2)+CONSTANTS.T./m.*u1;
vthetadot = -vr.*vtheta./r+CONSTANTS.T./m.*u2;
mdot   = -CONSTANTS.mdot*ones(size(t));
path = u(:,1).^2 + u(:,2).^2;

dae = [rdot vrdot vthetadot mdot path];
