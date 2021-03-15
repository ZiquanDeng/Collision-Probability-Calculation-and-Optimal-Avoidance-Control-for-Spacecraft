
function [dae] = spacecraftDae(sol);
global CONSTANTS;
t = sol.time;
x = sol.state;
u = sol.control;
p = sol.parameter;

% t0 = sol.initial.time;
% x0 = sol.initial.state;
% tf = sol.terminal.time;
% xf = sol.terminal.state;

% xf = sol.terminal.state;
% rf = xf(:,1:3);
% vf = xf(:,4:6);
% t0 = sol.initial.time;
% x0 = sol.initial.state;
% tf = sol.terminal.time;
ux = u(:,1);
uy = u(:,2);
uz = u(:,3);
r = x(:,1:3);
v = x(:,4:6);
g=9.8*10^-3;
Ft=0.2*g;
mu=3.986005e+5;
r1=6378.004+600;
w=sqrt(mu/r1^3);
rdot(:,1) = v(:,1);
rdot(:,2) = v(:,2);
rdot(:,3) = v(:,3);
vdot(:,1) = 2*w*v(:,3)+Ft*ux;
vdot(:,2) = -w^2*r(:,2)+Ft*uy;
vdot(:,3) = -2*w*v(:,1)+3*w^2*r(:,3)+Ft*uz;
path = (u(:,1).^2+u(:,2).^2+u(:,3).^2).^0.5;
dae = [rdot(:,1) rdot(:,2) rdot(:,3) vdot(:,1) vdot(:,2) vdot(:,3) path];





