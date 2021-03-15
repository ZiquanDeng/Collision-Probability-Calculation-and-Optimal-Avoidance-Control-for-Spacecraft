% ----------------------------------------------
% Minimum Time-to-Climb of a Supersonic Aircraft
% ----------------------------------------------
% This example is taken verbatim from the following reference:
% Betts, J. T., Practical Methods for Optimal Control Using
% Nonlinear Programming, SIAM Press, Philadelphia, 2001.
clear setup limits guess

% Initialize all of the data for the problem
load brysonMinimumClimbAeroData.mat;
global CONSTANTS;

% U.S. 1976 Standard Atmosphere Table
% Reference:  U.S. 1976 Standard Atmosphere, National Oceanographic
% and Atmospheric Administration, 1976.
% Column 1:  Altitude (m)
% Column 2:  Atmospheric Density (kg/m^3)
% Column 3:  Speed of Sound (m/s)
us1976 = [-2000    1.478e+00     3.479e+02
    0    1.225e+00     3.403e+02
    2000    1.007e+00     3.325e+02
    4000    8.193e-01     3.246e+02
    6000    6.601e-01     3.165e+02
    8000    5.258e-01     3.081e+02
    10000    4.135e-01     2.995e+02
    12000    3.119e-01     2.951e+02
    14000    2.279e-01     2.951e+02
    16000    1.665e-01     2.951e+02
    18000    1.216e-01     2.951e+02
    20000    8.891e-02     2.951e+02
    22000    6.451e-02     2.964e+02
    24000    4.694e-02     2.977e+02
    26000    3.426e-02     2.991e+02
    28000    2.508e-02     3.004e+02
    30000    1.841e-02     3.017e+02
    32000    1.355e-02     3.030e+02
    34000    9.887e-03     3.065e+02
    36000    7.257e-03     3.101e+02
    38000    5.366e-03     3.137e+02
    40000    3.995e-03     3.172e+02
    42000    2.995e-03     3.207e+02
    44000    2.259e-03     3.241e+02
    46000    1.714e-03     3.275e+02
    48000    1.317e-03     3.298e+02
    50000    1.027e-03     3.298e+02
    52000    8.055e-04     3.288e+02
    54000    6.389e-04     3.254e+02
    56000    5.044e-04     3.220e+02
    58000    3.962e-04     3.186e+02
    60000    3.096e-04     3.151e+02
    62000    2.407e-04     3.115e+02
    64000    1.860e-04     3.080e+02
    66000    1.429e-04     3.044e+02
    68000    1.091e-04     3.007e+02
    70000    8.281e-05     2.971e+02
    72000    6.236e-05     2.934e+02
    74000    4.637e-05     2.907e+02
    76000    3.430e-05     2.880e+02
    78000    2.523e-05     2.853e+02
    80000    1.845e-05     2.825e+02
    82000    1.341e-05     2.797e+02
    84000    9.690e-06     2.769e+02
    86000    6.955e-06     2.741e+02];

% Mtab is a table of Mach number values
Mtab   = [0; 0.2; 0.4; 0.6; 0.8; 1; 1.2; 1.4; 1.6; 1.8];
% alttab is a table of altitude values (in ft)
alttab = [0 5000 10000 15000 20000 25000 30000 40000 50000 70000];
% Convert altitude table to meters
alttab = 0.3048*alttab;
% Ttab is a table of aircraft thrust values (in lbf)
% Ttab is taken from Bryson's 1969 Journal of Aircraft paper (also
% Betts' book), but the table has been extended via linear
% extrapolation to fill in the "missing" data points.
Ttab = 1000*[24.2 24.0  20.3 17.3 14.5 12.2 10.2 5.7 3.4 0.1;
    28.0 24.6 21.1 18.1 15.2 12.8 10.7 6.5 3.9 0.2;
    28.3 25.2 21.9 18.7 15.9 13.4 11.2 7.3 4.4 0.4;
    30.8 27.2 23.8 20.5 17.3 14.7 12.3 8.1 4.9 0.8;
    34.5 30.3 26.6 23.2 19.8 16.8 14.1 9.4 5.6 1.1;
    37.9 34.3 30.4 26.8 23.3 19.8 16.8 11.2 6.8 1.4;
    36.1 38.0 34.9 31.3 27.3 23.6 20.1 13.4 8.3 1.7;
    36.1 36.6 38.5 36.1 31.6 28.1 24.2 16.2 10.0 2.2;
    36.1 35.2 42.1 38.7 35.7 32.0 28.1 19.3 11.9 2.9;
    36.1 33.8 45.7 41.3 39.8 34.6 31.1 21.7 13.3 3.1];
% Convert Thrust to Newtons
Ttab = 4.448222*Ttab;

% M2 is the Mach number used to compute the aerodynamic coefficients
M2         = [0 0.4 0.8 0.9 1.0 1.2 1.4 1.6 1.8];
Clalphatab = [3.44 3.44 3.44 3.58 4.44 3.44 3.01 2.86 2.44];
CD0tab     = [0.013 0.013 0.013 0.014 0.031 0.041 0.039 0.036 0.035];
etatab     = [0.54 0.54 0.54 0.75 0.79 0.78 0.89 0.93 0.93];

CONSTANTS.CDdat     = CDdat;
CONSTANTS.CLdat     = CLdat;
CONSTANTS.etadat    = etadat;
CONSTANTS.M         = Mtab;
CONSTANTS.M2        = M2;
CONSTANTS.alt       = alttab;
CONSTANTS.T         = Ttab;
CONSTANTS.Clalpha   = Clalphatab;
CONSTANTS.CD0       = CD0tab;
CONSTANTS.eta       = etatab;
CONSTANTS.ppCLalpha = polyfit(CONSTANTS.M2,CONSTANTS.Clalpha,8);
CONSTANTS.ppCD0     = polyfit(CONSTANTS.M2,CONSTANTS.CD0,8);
CONSTANTS.ppeta     = polyfit(CONSTANTS.M2,CONSTANTS.eta,8);
CONSTANTS.Re        = 6378145;
CONSTANTS.mu        = 3.986e14;
CONSTANTS.S         = 49.2386;
CONSTANTS.g0        = 9.80665;
CONSTANTS.Isp       = 1600;
CONSTANTS.H         = 7254.24;
CONSTANTS.rho0      = 1.225;
CONSTANTS.us1976    = us1976;

mass0 = 19050.864;
% Compute Scale Factors
if 0,
  scales.length = CONSTANTS.Re;
  scales.area = scales.length*scales.length;
  scales.volume = scales.area*scales.length;
  scales.speed = sqrt(CONSTANTS.mu/scales.length);
  scales.time = scales.length/scales.speed;
  scales.acceleration = scales.speed/scales.time;
  scales.mass = mass0;
  scales.force = scales.mass*scales.acceleration;
  scales.density = scales.mass/scales.volume;
  scales.gravparameter = scales.volume/scales.time^2;
else
  scales.length = 1;
  scales.area = 1;
  scales.volume = 1;
  scales.speed = 1;
  scales.time = 1;
  scales.acceleration = 1;
  scales.mass = 1;
  scales.force = 1;
  scales.density = 1;
  scales.gravparameter = 1;
end;

% Scale all quantities in problem
CONSTANTS.us1976(:,1) = us1976(:,1)/scales.length;
CONSTANTS.us1976(:,2) = us1976(:,2)/scales.density;
CONSTANTS.us1976(:,3) = us1976(:,3)/scales.speed;
CONSTANTS.alt         = CONSTANTS.alt/scales.length;
CONSTANTS.T           = CONSTANTS.T/scales.force;
CONSTANTS.Re          = CONSTANTS.Re/scales.length;
CONSTANTS.mu          = CONSTANTS.mu/scales.gravparameter;
CONSTANTS.S           = CONSTANTS.S/scales.area;
CONSTANTS.g0          = CONSTANTS.g0/scales.acceleration;
CONSTANTS.Isp         = CONSTANTS.Isp/scales.time;
CONSTANTS.H           = CONSTANTS.H/scales.length;
CONSTANTS.rho0        = CONSTANTS.rho0/scales.density;
[aa,mm]               = meshgrid(alttab,Mtab);
CONSTANTS.aa          = aa;
CONSTANTS.mm          = mm;

% Boundary conditions (taken from Betts 2001 & converted to SI units)
t0   = 0/scales.time;          % Initial time [scales.time]
alt0 = 0/scales.time;          % Initial altitude [scales.length]
altf = 19994.88/scales.length; % Final altitude [scales.length]
speed0 = 129.314/scales.speed; % Initial speed [scales.speed]
speedf = 295.092/scales.speed; % Final speed [scales.speed]
fpa0   = 0;                    % Initial flight path angle [rad]
fpaf   = 0;                    % Final flight path angle [rad]
mass0  = mass0/scales.mass;    % Initial mass [scales.mass]

% Bounds on variables (taken from Betts, 2001)
tfmin  = 100/scales.time;
tfmax  = 800/scales.time;
altmin = 0/scales.length;
altmax = 21031.2/scales.length;
speedmin = 5/scales.speed;
speedmax = 1000/scales.speed;
fpamin   = -40*pi/180;
fpamax   =  40*pi/180;
massmin  = 22/scales.mass;
massmax  = 20410/scales.mass;
% alphamin = -pi/2;
% alphamax = pi/2;
alphamin = -pi/4;
alphamax = pi/4;

iphase = 1;
% Bounds on initial and terminal values of time
% limits(iphase).meshPoints = [-1 1];
% limits(iphase).nodesPerInterval = [20];
limits(iphase).time.min = [t0 tfmin];
limits(iphase).time.max = [t0 tfmax];
limits(iphase).state.min(1,:) = [alt0 altmin altf];
limits(iphase).state.max(1,:) = [alt0 altmax altf];
limits(iphase).state.min(2,:) = [speed0 speedmin speedf];
limits(iphase).state.max(2,:) = [speed0 speedmax speedf];
limits(iphase).state.min(3,:) = [fpa0 fpamin fpaf];
limits(iphase).state.max(3,:) = [fpa0 fpamax fpaf];
limits(iphase).state.min(4,:) = [mass0 massmin massmin];
limits(iphase).state.max(4,:) = [mass0 massmax massmax];
limits(iphase).control.min = alphamin;
limits(iphase).control.max = alphamax;
limits(iphase).parameter.min = [];
limits(iphase).parameter.max = [];
guess(iphase).time = [0; 100];
guess(iphase).state(:,1) = [alt0; altf];
guess(iphase).state(:,2) = [speed0; speedf];
guess(iphase).state(:,3) = [fpa0; fpaf];
guess(iphase).state(:,4) = [mass0; mass0];
guess(iphase).control = [20; -20]*pi/180;
guess(iphase).parameter = [];

setup.name  = 'Bryson-Minimum-Time-to-Climb-Problem';
setup.funcs.cost = 'brysonMinimumClimbCost';
setup.funcs.dae = 'brysonMinimumClimbDae';
setup.funcs.link = '';
setup.limits = limits;
setup.guess = guess;
%================================================%
% WARNING:  AT THIS TIME THIS PROBLEM CAN ONLY   %
% BE SOLVED USING NUMERICAL DIFFERENTIATION!     %
% DO NOT SET "SETUP.DERIVATIVES" TO ANYTHING BUT %
% "finite-difference".                           %
%================================================%
setup.derivatives = 'finite-difference';
setup.autoscale = 'on';
% setup.tolerances = [1e-3 2e-3];
setup.mesh.tolerance = 1e-4;
setup.mesh.iteration = 10;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;

output = gpops(setup);

solution = output.solution;
solutionPlot = output.solutionPlot;
plotfigures;

%-----------------------------------------%
% End Function:  brysonMinimumClimbMain.m %
%-----------------------------------------%
