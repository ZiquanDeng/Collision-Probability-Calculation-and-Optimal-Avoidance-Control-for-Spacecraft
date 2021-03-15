%---------------------------------------------------%
% Classical Brachistochrone Problem:                %
%---------------------------------------------------%
% The problem solved here is given as follows:      %
%   Minimize t_f                                    %
% subject to the dynamic constraints                %
%    dx/dt =  v*sin(u)                              %
%    dy/dt = -v*cos(u)                              %
%    dv/dt =  g*cos(u)                              %
% and the boundary conditions                       %
%    x(0) = 0, y(0) = 0, v(0) = 0                   %
%    x(t_f) = 2, y(t_f) = -2, v(t_f) = FREE         %
%---------------------------------------------------%

clear all
close all
clc

global CONSTANTS 

CONSTANTS.g = 10;
x0 = 0;
y0 = 0;
v0 = 0;
xf = 2;
yf = -2;
xmin = -50;
xmax =  50;
ymin = -50;
ymax =   0;
vmin = xmin;
vmax = xmax;
param_min = [];
param_max = [];
path_min = [];
path_max = [];
event_min = [];
event_max = [];
duration_min = [];
duration_max = [];

iphase = 1;
limits(iphase).time.min         = [0 0];
limits(iphase).time.max         = [0 100];
limits(iphase).state.min(1,:)   = [x0 xmin xf];
limits(iphase).state.max(1,:)   = [x0 xmax xf];
limits(iphase).state.min(2,:)   = [y0 ymin yf];
limits(iphase).state.max(2,:)   = [y0 ymax yf];
limits(iphase).state.min(3,:)   = [v0 vmin vmin];
limits(iphase).state.max(3,:)   = [v0 vmax vmax];
limits(iphase).control.min(1,:) = -pi/2;
limits(iphase).control.max(1,:) =  pi/2;
limits(iphase).parameter.min    = param_min;
limits(iphase).parameter.max    = param_max;
limits(iphase).path.min         = path_min;
limits(iphase).path.max         = path_max;
limits(iphase).event.min        = event_min;
limits(iphase).event.max        = event_max;
limits(iphase).duration.min     = duration_min;
limits(iphase).duration.max     = duration_max;
guess(iphase).time              = [0; 20];
guess(iphase).state(:,1)        = [x0; xf];
guess(iphase).state(:,2)        = [y0; yf];
guess(iphase).state(:,3)        = [v0; 6];
guess(iphase).control(:,1)      = [0; 0];
% guess(iphase).control(:,2)      = [1; 1];
guess(iphase).parameter         = [];

setup.name  = 'Brachistochrone-Problem';
setup.funcs.cost = 'brachistochroneCost';
setup.funcs.dae = 'brachistochroneDae';
setup.limits = limits;
setup.guess = guess;
setup.linkages = [];
setup.derivatives = 'finite-difference';
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-6;
setup.mesh.iteration = 20;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;

[output,gpopsHistory] = gpops(setup);
solution = output.solution;
solutionPlot = output.solutionPlot;
