%--------------------------------------------------------%
% Classical Hohmann Transfer Between Two Circular Orbits %
%--------------------------------------------------------%
clear all
clc

global CONSTANTS
CONSTANTS.mu = 3.986e14;
CONSTANTS.R_e = 6378145;

% Define Boundary Conditions
h0 = 200000;  % Initial Altitude  = 200 km
hf = 400000;  % Terminal Altitude = 400 km
r0 = CONSTANTS.R_e + h0;
a0 = r0;
e0 = 0;
omega0 = 0;
cap_omega0 = 0;
i0 = 0;
f0 = 0;
elmts = [a0,e0,i0,omega0,cap_omega0,f0];
[r_eci_0,v_eci_0minus] = hohmannTransferOe2rv(elmts,CONSTANTS.mu);
af = hf + CONSTANTS.R_e;
ef = 0;
CONSTANTS.v0 = v_eci_0minus;
Rf = hf + CONSTANTS.R_e;
vf = sqrt(CONSTANTS.mu/Rf);
gf = 0;
t0 = 0;
tfmin = 1000;
tfmax = 3000;
xmin = -1.2*(hf + CONSTANTS.R_e);
xmax = -xmin;
ymin = -10*(hf + CONSTANTS.R_e);
ymax = -ymin;
zmin = -10*(hf + CONSTANTS.R_e);
zmax = -zmin;
vxmin = -10*sqrt(CONSTANTS.mu/r0);
vxmax = -vxmin;
vymin = vxmin;
vymax = vxmax;
vzmin = 0;
vzmax = 0;
x0 = r_eci_0(1);
y0 = r_eci_0(2);
z0 = r_eci_0(3);
vx0 = CONSTANTS.v0(1);
vy0 = CONSTANTS.v0(2);
vz0 = CONSTANTS.v0(3);
numin = 50*pi/180;
numax = 2*pi;
nuf = pi;

iphase = 1;
limits(iphase).time.min = [t0 tfmin];
limits(iphase).time.max = [t0 tfmax];
limits(iphase).state.min(1,:)   = [x0 xmin xmin];
limits(iphase).state.max(1,:)   = [x0 xmax xmax];
limits(iphase).state.min(2,:)   = [y0 ymin ymin];
limits(iphase).state.max(2,:)   = [y0 ymax ymax];
limits(iphase).state.min(3,:)   = [z0 zmin zmin];
limits(iphase).state.max(3,:)   = [z0 zmax zmax];
limits(iphase).state.min(4,:)   = [vxmin vxmin vxmin];
limits(iphase).state.max(4,:)   = [vxmax vxmax vxmax];
limits(iphase).state.min(5,:)   = [vymin vymin vymin];
limits(iphase).state.max(5,:)   = [vymax vymax vymax];
limits(iphase).state.min(6,:)   = [vzmin vzmin vzmin];
limits(iphase).state.max(6,:)   = [vzmax vzmax vzmax];
limits(iphase).control.min = [];
limits(iphase).control.max = [];
limits(iphase).parameter.min    = -1000*ones(6,1);
limits(iphase).parameter.max    = 1000*ones(6,1);

limits(iphase).path.min    = [];
limits(iphase).path.max    = [];
limits(iphase).event.min    = [];
limits(iphase).event.max    = [];
limits(iphase).event.min   = [0; 0; 0; Rf; vf; gf];
limits(iphase).event.max   = [0; 0; 0; Rf; vf; gf];

tfGuess = 2500;
guess(iphase).time =  [t0; tfGuess];
guess(iphase).state(:,1) = [x0; x0];
guess(iphase).state(:,2) = [y0; y0];
guess(iphase).state(:,3) = [z0; z0];
guess(iphase).state(:,4) = [vx0; vx0];
guess(iphase).state(:,5) = [vy0; vy0];
guess(iphase).state(:,6) = [vz0; vz0];
guess(iphase).parameter = [100; 100; 100; 100; 100; 100];

setup.autoscale = 'on';
setup.name = 'Hohmann-Transfer';
setup.funcs.cost = 'hohmannTransferCost';
setup.funcs.dae = 'hohmannTransferDae';
setup.funcs.event = 'hohmannTransferEvent';
setup.derivatives = 'automatic-intlab';
setup.mesh.tolerance = 1e-6;
setup.mesh.iteration = 20;
setup.limits = limits;
setup.linkages = [];
setup.guess = guess;
output = gpops(setup);
solution = output.solution;

