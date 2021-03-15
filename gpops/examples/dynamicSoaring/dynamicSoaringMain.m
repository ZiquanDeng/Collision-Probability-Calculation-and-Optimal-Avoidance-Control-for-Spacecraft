%-----------------------------------------------------------------------%
% BEGIN script dynamicSoaringMain.m                                     %
%-----------------------------------------------------------------------%
% This example is taken from the following reference:                   %
% Zhao, Y. J., "Optimal Pattern of Glider Dynamic Soaring," Optimal     %
% Control Applications and Methods, Vol. 25, 2004, pp. 67-89            %
%                                                                       %
% The optimal control problem is described as follows.                  %
% Minimize beta                                                         %
% subject to the dynamic constraints                                    %
%   xdot = V*cos(gamma)*sin(psi)+Wx                                     %
%   ydot = V*cos(gamma)*cos(psi)                                        %
%   hdot = V*sin(gamma)                                                 %
%   m*vdot=-D-mg*sin(gamma)-m*Wxdot*cos(gamma)*sin(psi)                 %
%   m*V*gammadot=L*cos(sigma)-mg*cos(gamma)+m*Wxdot*sin(gamma)*sin(psi) %
%   m*V*cos(gamma)*psidot=L*sin(sigma)-m*Wxdot*cos(psi)                 %
% the inequality path constraint                                        %
%   -2 <= L/(mg) <= 5                                                   %
% and the boundary conditions                                           %
%   x(0)=0, y(0)=0, z(0)=0, x(tf)=x(0), y(tf)=y(0), h(tf)=h(0),         %
%   v(tf)=v(0), gamma(tf)=gamma(0), psi(tf) = psi(0)-2*pi               %
%-----------------------------------------------------------------------%
clear all
clc

global CONSTANTS
CONSTANTS.rho=0.002378;
CONSTANTS.CD0 = 0.00873; 
CONSTANTS.K= 0.045;
CONSTANTS.g=32.2;
CONSTANTS.m=5.6;
CONSTANTS.S=45.09703;
CONSTANTS.mu = 3.986e14;

CONSTANTS.mgos=CONSTANTS.m*CONSTANTS.g/CONSTANTS.S;
CONSTANTS.Emax=(1/(4*CONSTANTS.K*CONSTANTS.CD0))^0.5;
CONSTANTS.W0=0;
CONSTANTS.lmin = -2;
CONSTANTS.lmax = 5; %load factor constraint

% CONSTANTS.ls=6378145;
% CONSTANTS.vs=sqrt(CONSTANTS.mu/CONSTANTS.ls);
% CONSTANTS.ts=CONSTANTS.ls/CONSTANTS.vs;

CONSTANTS.ls = 1;
CONSTANTS.ts = 1;
CONSTANTS.vs = 1;

x0 = 0;
y0 = 0;
z0 = 0;
r0 = 0;
rf = 0;
v0=100;
v0=100/CONSTANTS.vs;

xmin = -1/CONSTANTS.ls*1e3;
xmax =  1/CONSTANTS.ls*1e3;
ymin = -1/CONSTANTS.ls*1e3;
ymax =  1/CONSTANTS.ls*1e3;
zmin =    1/CONSTANTS.ls*0;
zmax =  1/CONSTANTS.ls*1e3;
vmin =   1/CONSTANTS.vs*10;
vmax =  1/CONSTANTS.vs*350;
rmin = -75/90*pi/2;
rmax =  75/90*pi/2;
psimin = -3*pi;
psimax = 0.5*pi;

betamin =0.005;
betamax = 0.15;

CLmin = -0.5;
CLmin =  0;
CLmax = 1.5;
Phimin = -75/180*pi;
Phimax =  75/180*pi;
t0 = 0;
tfmin = 1/CONSTANTS.ts;
tfmax = 30/CONSTANTS.ts;

% Phase 1 Information
iphase = 1;
limits(iphase).time.min = [t0 tfmin];
limits(iphase).time.max = [t0 tfmax];
limits(iphase).state.min(1,:) = [x0 xmin x0];
limits(iphase).state.max(1,:) = [x0 xmax x0];
limits(iphase).state.min(2,:) = [y0 ymin y0];
limits(iphase).state.max(2,:) = [y0 ymax y0];
limits(iphase).state.min(3,:) = [z0 zmin z0];
limits(iphase).state.max(3,:) = [z0 zmax z0];
limits(iphase).state.min(4,:) = [vmin vmin vmin];
limits(iphase).state.max(4,:) = [vmax vmax vmax];
limits(iphase).state.min(5,:) = [rmin rmin rmin];
limits(iphase).state.max(5,:) = [rmax rmax rmax];
limits(iphase).state.min(6,:) = [psimin psimin psimin];
limits(iphase).state.max(6,:) = [psimax psimax psimax];

if 1,
limits(iphase).control.min(1,:) = CLmin;
limits(iphase).control.max(1,:) = CLmax;
limits(iphase).control.min(2,:) = Phimin;
limits(iphase).control.max(2,:) = Phimax;
else
limits(iphase).control.min(1,:) = -CLmax;
limits(iphase).control.max(1,:) =  CLmax;
limits(iphase).control.min(2,:) = -CLmax;
limits(iphase).control.max(2,:) =  CLmax;
end;

limits(iphase).parameter.min(1,1) = [betamin];
limits(iphase).parameter.max(1,1) = [betamax];

if 1
  limits(iphase).path.min = [CONSTANTS.lmin];
  limits(iphase).path.max = [CONSTANTS.lmax];
else
  limits(iphase).path.min = [CONSTANTS.lmin; 0];
  limits(iphase).path.max = [CONSTANTS.lmax; CLmax^2];
end;
limits(iphase).event.min     = [0 ; 0 ; 0];
limits(iphase).event.max     = [0 ; 0 ; 0];

limits(iphase).duration.min = [];
limits(iphase).duration.max = [];

% N=100;
% CL0=CLmax;
% basetime=linspace(0,1,N)';
% xguess=zeros(N,1);
% yguess=zeros(N,1);
% zguess=zeros(N,1);
% vguess=v0*(2+cos(2*pi*basetime));
% rguess=pi/6*sin(2*pi*basetime);
% psiguess=-2*pi*sin(1/2*pi*basetime);
% CLguess=CL0*ones(N,1);
% phiguess= zeros(N,1);
% u1guess = zeros(N,1);
% u2guess = zeros(N,1);
% betaguess=0.08;
% loadguess=0.2;

N=100;
CL0=CLmax;
basetime=linspace(0,24,N)';
xguess=500*cos(2*pi*basetime/24)-500;
yguess=300*sin(2*pi*basetime/24);
zguess=-400*cos(2*pi*basetime/24)+400;
vguess=0.8*v0*(1.5+cos(2*pi*basetime/24));
rguess=pi/6*sin(2*pi*basetime/24);
psiguess=-1-basetime/4;
CLguess=CL0*ones(N,1)/3;
phiguess= -ones(N,1);
%u1guess = zeros(N,1);
%u2guess = zeros(N,1);
betaguess=0.08;
%loadguess=0.2;

guess(iphase).time            = [basetime];
guess(iphase).state(:,1)      = [xguess];
guess(iphase).state(:,2)      = [yguess];
guess(iphase).state(:,3)      = [zguess];
guess(iphase).state(:,4)      = [vguess];
guess(iphase).state(:,5)      = [rguess];
guess(iphase).state(:,6)      = [psiguess];
guess(iphase).control(:,1)    = [CLguess];
guess(iphase).control(:,2)    = [phiguess];
guess(iphase).parameter(1,1)  = [betaguess];


connections = [];
setup.name  = 'DynamicSoaring-Problem';
setup.funcs.cost = 'dynamicSoaringCost';
setup.funcs.dae = 'dynamicSoaringDae';
setup.funcs.event = 'dynamicSoaringEvent';
setup.limits = limits;
setup.guess = guess;
setup.linkages = [];
setup.derivatives = 'finite-difference';
setup.autoscale = 'on';   %If you need autoscaling turn on [on, off].
setup.mesh.tolerance = 1e-6;
setup.mesh.iteration = 20;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;

output = gpops(setup);
solution = output.solution;
solutionPlot = output.solutionPlot;
