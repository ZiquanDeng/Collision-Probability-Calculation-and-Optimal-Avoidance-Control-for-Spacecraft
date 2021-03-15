%-----------------------------------------------------%
% Goddard Rocket Problem                              %
%-----------------------------------------------------%
% This problem is taken from the following reference: %
% Bryson, A. E. and Ho, Y-C, Applied Optiml Control,  %
% Hemisphere Publishing, 1975.                        %
%-----------------------------------------------------%
clear all
clc

global CONSTANTS

h0 = 0;
v0 = 0;
m0 = 3;
mf = 1;

hmin = 0;
hmax = 30000;
vmin = -1500;
vmax = 1500;
mmin = 0.2*m0;
mmax = m0;
t0 = 0;
tfMin = 0;
tfMax = 500;

CONSTANTS.g0 = 32.174;
CONSTANTS.rho0 = 0.002378;
CONSTANTS.H = 23800;
CONSTANTS.csqrd = 3.264*CONSTANTS.g0*CONSTANTS.H;
CONSTANTS.c = sqrt(CONSTANTS.csqrd);
Tmax = 2*m0*CONSTANTS.g0;
CONSTANTS.dragk = 0.7110*Tmax/CONSTANTS.csqrd;

iphase = 1;
limits(iphase).time.min = [t0 tfMin];
limits(iphase).time.max = [t0 tfMax];
limits(iphase).state.min(1,:)   = [h0 hmin hmin];
limits(iphase).state.max(1,:)   = [h0 hmax hmax];
limits(iphase).state.min(2,:)   = [v0 vmin vmin];
limits(iphase).state.max(2,:)   = [v0 vmax vmax];
limits(iphase).state.min(3,:)   = [m0 mmin mmin];
limits(iphase).state.max(3,:)   = [m0 mmax mmax];
limits(iphase).control.min(1,:) = Tmax;
limits(iphase).control.max(1,:) = Tmax;
limits(iphase).parameter.min    = [];
limits(iphase).parameter.max    = [];
limits(iphase).path.min    = [];
limits(iphase).path.max    = [];
limits(iphase).event.min   = [];
limits(iphase).event.max   = [];
limits(iphase).duration.min = 1;
limits(iphase).duration.max = 1000;
guess(iphase).time =  [tfMin; tfMax];
guess(iphase).state(:,1) = [h0; h0];
guess(iphase).state(:,2) = [v0; v0];
guess(iphase).state(:,3) = [m0; mf];
guess(iphase).control(:,1) = [0; Tmax];
guess(iphase).parameter = []; 

iphase = 2;
limits(iphase).time.min = [tfMin tfMin];
limits(iphase).time.max = [tfMax tfMax];
limits(iphase).state.min(1,:)   = [hmin hmin hmin];
limits(iphase).state.max(1,:)   = [hmax hmax hmax];
limits(iphase).state.min(2,:)   = [vmin vmin vmin];
limits(iphase).state.max(2,:)   = [vmax vmax vmax];
limits(iphase).state.min(3,:)   = [mmin mmin mmin];
limits(iphase).state.max(3,:)   = [mmax mmax mmax];
limits(iphase).control.min(1,:) = 0;
limits(iphase).control.max(1,:) = Tmax;
limits(iphase).parameter.min    = [];
limits(iphase).parameter.max    = [];
limits(iphase).path.min    = [];
limits(iphase).path.max    = [];
limits(iphase).path.min    = 0;
limits(iphase).path.max    = 0;
limits(iphase).event.min   = [0];
limits(iphase).event.max   = [0];
limits(iphase).duration.min = 1;
limits(iphase).duration.max = 1000;
guess(iphase).time =  [tfMin; tfMax];
guess(iphase).state(:,1) = [h0; hmax];
guess(iphase).state(:,2) = [vmax; vmax];
guess(iphase).state(:,3) = [m0; mf];
guess(iphase).control(:,1) = [0; Tmax];
guess(iphase).parameter = [];

iphase = 3;
limits(iphase).time.min = [tfMin tfMin];
limits(iphase).time.max = [tfMax tfMax];
limits(iphase).state.min(1,:)   = [hmin hmin hmin];
limits(iphase).state.max(1,:)   = [hmax hmax hmax];
limits(iphase).state.min(2,:)   = [vmin vmin vmin];
limits(iphase).state.max(2,:)   = [vmax vmax vmax];
limits(iphase).state.min(3,:)   = [mmin mmin mf];
limits(iphase).state.max(3,:)   = [mmax mmax mf];
limits(iphase).control.min(1,:) = 0;
limits(iphase).control.max(1,:) = 0;
limits(iphase).parameter.min    = [];
limits(iphase).parameter.max    = [];
limits(iphase).path.min    = [];
limits(iphase).path.max    = [];
limits(iphase).event.min   = [];
limits(iphase).event.max   = [];
limits(iphase).duration.min = 1;
limits(iphase).duration.max = 1000;
guess(iphase).time =  [tfMin; tfMax];
guess(iphase).state(:,1) = [hmax; hmax];
guess(iphase).state(:,2) = [vmax; vmax];
guess(iphase).state(:,3) = [m0; mf];
guess(iphase).control(:,1) = [0; Tmax];
guess(iphase).parameter = [];

ipair = 1;
linkages(ipair).left.phase = 1;
linkages(ipair).right.phase = 2;
linkages(ipair).min = [0; 0; 0];
linkages(ipair).max = [0; 0; 0];

ipair = 2;
linkages(ipair).left.phase = 2;
linkages(ipair).right.phase = 3;
linkages(ipair).min = [0; 0; 0];
linkages(ipair).max = [0; 0; 0];

setup.name = 'Goddard-Rocket';
setup.funcs.cost = 'goddardRocketCost';
setup.funcs.dae = 'goddardRocketDae';
setup.funcs.event = 'goddardRocketEvent';
setup.funcs.link = 'goddardRocketLink';
setup.limits = limits;
setup.derivatives = 'finite-difference';
setup.guess = guess;
setup.linkages = linkages;
setup.autoscale = 'off';
setup.mesh.tolerance = 1e-6;
setup.mesh.iteration = 20;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;

[output,gpopsHistory] = gpops(setup);
solution = output.solution;
solutionPlot = output.solutionPlot;

%--------------------------------%
% End File:  goddardRocketMain.m %
%--------------------------------%
