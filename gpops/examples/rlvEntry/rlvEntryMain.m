% --------------------------------------
% Reusable Launch Vehicle Entry Example
% --------------------------------------
% This example is taken verbatim from the following reference:
%   Betts, J. T., Practical Methods for Optimal Control Using
%   Nonlinear Programming, SIAM Press, Philadelphia, 2001.

clear all
close all
clc

global constants

constants.Re = 20902900;         % Equatorial Radius of Earth (ft)
constants.S  = 2690;             % Vehicle Reference Area (ft^2)
constants.cl(1) = -0.2070;       % Parameters for lift coefficient
constants.cl(2) = 1.6756;        % cl(1) and cl(2)
constants.cd(1) = 0.0785;        % Parameters for drag coefficient
constants.cd(2) = -0.3529;       % cd(1), cd(2), and cd(3)
constants.cd(3) = 2.0400;
constants.b(1)  = 0.07854;       % Coefficient for heating rate
constants.b(2)  = -0.061592;     % b(1),  b(2), and b(3)
constants.b(3)  = 0.00621408;
constants.H     = 23800;         % Density Scale Height (ft)
constants.al(1) = -0.20704;    
constants.al(2) = 0.029244;
constants.rho0  = 0.002378;      % Sea Level Atmospheric Density (slug/ft^3)
constants.mu    = 1.4076539e16;  % Earth Gravitational Parameter (ft^^3/s^2) 
constants.mass  = 6309.433;      % 

scales.length = constants.Re;
scales.speed = sqrt(constants.mu/constants.Re);
scales.time = scales.length/scales.speed;
scales.acceleration = scales.speed/scales.time;
scales.area = scales.length*scales.length;
scales.volume = scales.area*scales.length;
scales.mass = constants.mass;
scales.density = scales.mass/scales.volume;
scales.mu = scales.volume/scales.time^2;

if 1,
scales.length = 1;
scales.speed = 1;
scales.time = 1;
scales.acceleration = 1;
scales.area = 1;
scales.volume = 1;
scales.mass = 1;
scales.density = 1;
scales.mu = 1;
end;

constants.Re = constants.Re/scales.length;   
constants.S  = constants.S/scales.area;       
constants.cl(1) = -0.2070; 
constants.cl(2) = 1.6756;  
constants.cd(1) = 0.0785;  
constants.cd(2) = -0.3529; 
constants.cd(3) = 2.0400;
constants.b(1)  = 0.07854; 
constants.b(2)  = -0.061592;
constants.b(3)  = 0.00621408;
constants.H     = constants.H/scales.length;     
constants.al(1) = -0.20704;    
constants.al(2) = 0.029244;
constants.rho0  = 0.002378/scales.density;     
constants.mu    = constants.mu/scales.mu;
constants.mass  = constants.mass/scales.mass;    


t0 = 0/scales.time;
alt0 = 260000/scales.length;
rad0 = alt0+constants.Re;
lon0 = 0;
lat0 = 0;
speed0 = 25600/scales.speed;
fpa0   = -1*pi/180;
azi0   = 90*pi/180;

tf = 50/scales.time;
altf = 80000/scales.length;
radf = altf+constants.Re;
lonf = 90*pi/180;
latf = 30*pi/180;
speedf = 2500/scales.speed;
fpaf   = -5*pi/180;
azif   = -90*pi/180;

t0min = 0/scales.time;
t0max = 0/scales.time;
tfmin = 0/scales.time;
tfmax = 3000/scales.time;
radmin = constants.Re;
radmax = rad0;
lonmin = -180*pi/180;
lonmax = -lonmin;
latmin = -70*pi/180;
latmax = -latmin;
speedmin = 10/scales.speed;
speedmax = 45000/scales.speed;
fpamin = -80*pi/180;
fpamax =  80*pi/180;
azimin = -180*pi/180;
azimax =  180*pi/180;
aoamin = -90*pi/180;
aoamax = -aoamin;
bankmin = -90*pi/180;
bankmax =   1*pi/180;

iphase = 1;

limits(iphase).time.min    = [t0min tfmin];
limits(iphase).time.max    = [t0max tfmax];
limits(iphase).state.min(1,:) = [rad0 radmin radf];
limits(iphase).state.max(1,:) = [rad0 radmax radf];
limits(iphase).state.min(2,:) = [lon0 lonmin lonmin];
limits(iphase).state.max(2,:) = [lon0 lonmax lonmax];
limits(iphase).state.min(3,:) = [lat0 latmin latmin];
limits(iphase).state.max(3,:) = [lat0 latmax latmax];
limits(iphase).state.min(4,:) = [speed0 speedmin speedf];
limits(iphase).state.max(4,:) = [speed0 speedmax speedf];
limits(iphase).state.min(5,:) = [fpa0 fpamin fpaf];
limits(iphase).state.max(5,:) = [fpa0 fpamax fpaf];
limits(iphase).state.min(6,:) = [azi0 azimin azimin];
limits(iphase).state.max(6,:) = [azi0 azimax azimax];
limits(iphase).control.min(1,:) = aoamin;
limits(iphase).control.max(1,:) = aoamax;
limits(iphase).control.min(2,:) = bankmin;
limits(iphase).control.max(2,:) = bankmax;
limits(iphase).parameter.min = [];
limits(iphase).parameter.max = [];
limits(iphase).path.min      = [];
limits(iphase).path.max      = [];
limits(iphase).event.min     = [];
limits(iphase).event.max     = [];
guess(iphase).time            = [0; 1000];
guess(iphase).state(:,1)     = [rad0; radf];
guess(iphase).state(:,2)     = [lon0; lon0];
guess(iphase).state(:,3)     = [lat0; lat0];
guess(iphase).state(:,4)     = [speed0; speedf];
guess(iphase).state(:,5)     = [fpa0; fpaf];
guess(iphase).state(:,6)     = [azi0; azif];
guess(iphase).control(:,1)   = [0; 0];
guess(iphase).control(:,2)   = [0; 0];
guess(iphase).parameter      = [];

% clear guess
% load guess.mat

setup.name = 'RLV-Entry';
setup.limits = limits;
setup.guess = guess;
setup.funcs.cost = 'rlvEntryCost';
setup.funcs.dae  = 'rlvEntryDae';
setup.linkages = [];
setup.derivatives = 'automatic-intlab';
setup.autoscale = 'on';
% setup.tolerances = [1e-6 2e-6];
setup.mesh.tolerance = 1e-6;
setup.mesh.iteration = 10;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 10;
[output,gpopsHistory] = gpops(setup);
solution = output.solution;
solutionPlot = output.solutionPlot;

%------------------------------%
% END: function rlvEntryMain.m %
%------------------------------%
