function dae = orbitRaisingMaxRadiusDae(soldae,iphase);

global CONSTANTS

t = soldae.time;
x = soldae.state;
u = soldae.control;
p = soldae.parameter;

r      = x(:,1);
theta  = x(:,2);
vr     = x(:,3);
vtheta = x(:,4);
m         = CONSTANTS.m0-CONSTANTS.mdot.*t;
a         = CONSTANTS.T./m;
rdot      = vr;
thetadot  = vtheta./r;
vrdot     = vtheta.^2./r-CONSTANTS.mu./r.^2+a.*u(:,1);
vthetadot = -vr.*vtheta./r+a.*u(:,2); 
path      = u(:,1).^2+u(:,2).^2;
dae       = [rdot thetadot vrdot vthetadot path];
