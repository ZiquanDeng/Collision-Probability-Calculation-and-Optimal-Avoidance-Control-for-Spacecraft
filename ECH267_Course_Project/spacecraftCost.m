
function [Mayer,Lagrange]=spacecraftCost(sol);
global CONSTANTS;
t0 = sol.initial.time;
% x0 = sol.initial.state;
tf = sol.terminal.time;
xf = sol.terminal.state;
t  = sol.time;
% x  = sol.state;
u = sol.control;
% u2 = sol.control;
% u3 = sol.control;
u = sol.control;
p  = sol.parameter;
g=9.8*10^-3;
Ft=0.3*g;
ux = u(:,1);
uy = u(:,2);
uz = u(:,3);
 Mayer = zeros(size(t0));
% Lagrange = [u(:,1) u(:,2) u(:,3)]';

% Mayer =ux.^2+uy.^2+uz.^2
% Mayer =[ux uy uz]*[ux uy uz].';
% Mayer =tf;
% Mayer = ux.'*ux+uy.'*uy+uz.'*uz;
% 
% Mayer =[u(1:3)]*[u(1:3)].';
% Mayer =ux'.*ux+uy'.*uy+uz'.*uz
% Lagrange = zeros(size(t));
% Lagrange = [];

% Mayer = tf;
% Lagrange = zeros(size(t));
Lagrange = ux.^2+uy.^2+uz.^2;


