function [ode path] = tumorAntiAngiogenesisDae(soldae)

global CONSTANTS

t = soldae.time;
x = soldae.state;
u = soldae.control;
p = soldae.parameter;

p = x(:,1);
q = x(:,2);
y = x(:,3);

pdot = -CONSTANTS.zeta.*p.*log(p./q);
qdot = q.*(CONSTANTS.b - (CONSTANTS.mu + (CONSTANTS.d*(p.^(2/3))) + CONSTANTS.G.*u));
ydot = u;

ode = [pdot qdot ydot];
path = [];

