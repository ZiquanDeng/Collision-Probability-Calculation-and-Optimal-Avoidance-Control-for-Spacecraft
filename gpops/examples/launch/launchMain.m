%----------------------------------------------------------------------------%
%                 Multiple-Stage Launch Vehicle Ascent Example               %
% ---------------------------------------------------------------------------%
% This example can be found in the following reference:                      %
%   Rao, A. V., Benson, D. A., Darby, C. L., Patterson, M. A., Francolin,    % 
%   C., Sanders, I., and Huntington, G. T., “Algorithm 902: GPOPS, A MATLAB  %
%   Software for Solving Multiple-Phase Optimal Control Problems Using The   %
%   Gauss Pseudospectral Method,” ACM Transactions on Mathematical Software, %
%   Vol. 37, No. 2, April-June 2010, Article No. 22, pp. 1-39.               %
% ---------------------------------------------------------------------------%
clear all
clc

earthRadius         = 6378145;
gravParam           = 3.986012e14;
initialMass         = 301454;
earthRotRate        = 7.29211585e-5;
seaLevelDensity     = 1.225;
densityScaleHeight  = 7200;
g0                  = 9.80665;

scales.length       = earthRadius;
scales.speed        = sqrt(gravParam/scales.length);
scales.time         = scales.length/scales.speed;
scales.acceleration = scales.speed/scales.time;
scales.mass         = initialMass;
scales.force        = scales.mass*scales.acceleration;
scales.area         = scales.length^2;
scales.volume       = scales.area.*scales.length;
scales.density      = scales.mass/scales.volume;
scales.gravparam    = scales.acceleration*scales.length^2;

omega               = earthRotRate*scales.time;
auxdata.omega_matrix = omega*[0 -1 0;1 0 0;0 0 0];
auxdata.mu          = gravParam/scales.gravparam;
auxdata.cd          = 0.5;
auxdata.sa          = 4*pi/scales.area;
auxdata.rho0        = seaLevelDensity/scales.density;
auxdata.H           = densityScaleHeight/scales.length;
auxdata.Re          = earthRadius/scales.length;
auxdata.g0          = g0/scales.acceleration;


global CONSTANTS

omega               = earthRotRate*scales.time;
CONSTANTS.omega_matrix = omega*[0 -1 0;1 0 0;0 0 0];
CONSTANTS.mu          = gravParam/scales.gravparam;
CONSTANTS.cd          = 0.5;
CONSTANTS.sa          = 4*pi/scales.area;
CONSTANTS.rho0        = seaLevelDensity/scales.density;
CONSTANTS.H           = densityScaleHeight/scales.length;
CONSTANTS.Re          = earthRadius/scales.length;
CONSTANTS.g0          = g0/scales.acceleration;

lat0 = 28.5*pi/180;               % Geocentric Latitude of Cape Canaveral
x0 = CONSTANTS.Re*cos(lat0);      % x component of initial position
z0 = CONSTANTS.Re*sin(lat0);      % z component of initial position
y0 = 0;
r0 = [x0; y0; z0];
v0 = CONSTANTS.omega_matrix*r0;

bt_srb = 75.2/scales.time;
bt_first = 261/scales.time;
bt_second = 700/scales.time;

t0 = 0/scales.time;
t1 = 75.2/scales.time;
t2 = 150.4/scales.time;
t3 = 261/scales.time;
t4 = 961/scales.time;

m_tot_srb     = 19290/scales.mass;
m_prop_srb    = 17010/scales.mass;
m_dry_srb     = m_tot_srb-m_prop_srb;
m_tot_first   = 104380/scales.mass;
m_prop_first  = 95550/scales.mass;
m_dry_first   = m_tot_first-m_prop_first;
m_tot_second  = 19300/scales.mass;
m_prop_second = 16820/scales.mass;
m_dry_second  = m_tot_second-m_prop_second;
m_payload     = 4164/scales.mass;
thrust_srb    = 628500/scales.force;
thrust_first  = 1083100/scales.force;
thrust_second = 110094/scales.force;
mdot_srb      = m_prop_srb/bt_srb;
ISP_srb       = thrust_srb/(CONSTANTS.g0*mdot_srb);
mdot_first    = m_prop_first/bt_first;
ISP_first     = thrust_first/(CONSTANTS.g0*mdot_first);
mdot_second   = m_prop_second/bt_second;
ISP_second     = thrust_second/(CONSTANTS.g0*mdot_second);

af = 24361140/scales.length;
ef = 0.7308;
incf = 28.5*pi/180;
Omf = 269.8*pi/180;
omf = 130.5*pi/180;
nuguess = 0;
cosincf = cos(incf);
cosOmf = cos(Omf);
cosomf = cos(omf);
oe = [af ef incf Omf omf nuguess];
[rout,vout] = launchoe2rv(oe,CONSTANTS.mu);
rout = rout';
vout = vout';

m10 = m_payload+m_tot_second+m_tot_first+9*m_tot_srb;
m1f = m10-(6*mdot_srb+mdot_first)*t1;
m20 = m1f-6*m_dry_srb;
m2f = m20-(3*mdot_srb+mdot_first)*(t2-t1);
m30 = m2f-3*m_dry_srb;
m3f = m30-mdot_first*(t3-t2);
m40 = m3f-m_dry_first;
m4f = m_payload;

CONSTANTS.thrust_srb    = thrust_srb;
CONSTANTS.thrust_first  = thrust_first;
CONSTANTS.thrust_second = thrust_second;
CONSTANTS.ISP_srb       = ISP_srb;
CONSTANTS.ISP_first     = ISP_first;
CONSTANTS.ISP_second    = ISP_second;

rmin = -2*CONSTANTS.Re;
rmax = -rmin;
vmin = -10000/scales.speed;
vmax = -vmin;

iphase = 1;
limits(iphase).time.min = [t0 t1];
limits(iphase).time.max = [t0 t1];
limits(iphase).state.min(1,:) = [r0(1) rmin rmin];
limits(iphase).state.max(1,:) = [r0(1) rmax rmax];
limits(iphase).state.min(2,:) = [r0(2) rmin rmin];
limits(iphase).state.max(2,:) = [r0(2) rmax rmax];
limits(iphase).state.min(3,:) = [r0(3) rmin rmin];
limits(iphase).state.max(3,:) = [r0(3) rmax rmax];
limits(iphase).state.min(4,:) = [v0(1) vmin vmin];
limits(iphase).state.max(4,:) = [v0(1) vmax vmax];
limits(iphase).state.min(5,:) = [v0(2) vmin vmin];
limits(iphase).state.max(5,:) = [v0(2) vmax vmax];
limits(iphase).state.min(6,:) = [v0(3) vmin vmin];
limits(iphase).state.max(6,:) = [v0(3) vmax vmax];
limits(iphase).state.min(7,:) = [m10 m1f m1f];
limits(iphase).state.max(7,:) = [m10 m10 m10];
limits(iphase).control.min(1,:) = -1;
limits(iphase).control.max(1,:) =  1;
limits(iphase).control.min(2,:) = -1;
limits(iphase).control.max(2,:) =  1;
limits(iphase).control.min(3,:) = -1;
limits(iphase).control.max(3,:) =  1;
limits(iphase).parameter.min    = [];
limits(iphase).parameter.max    = [];
limits(iphase).path.min    = 1;
limits(iphase).path.max    = 1;
guess(iphase).time = [t0; t1];
guess(iphase).state(:,1) = [r0(1); r0(1)];
guess(iphase).state(:,2) = [r0(2); r0(2)];
guess(iphase).state(:,3) = [r0(3); r0(3)];
guess(iphase).state(:,4) = [v0(1); v0(1)];
guess(iphase).state(:,5) = [v0(2); v0(2)];
guess(iphase).state(:,6) = [v0(3); v0(3)];
guess(iphase).state(:,7) = [m10; m1f];
guess(iphase).control(:,1) = [0; 0];
guess(iphase).control(:,2) = [1; 1];
guess(iphase).control(:,3) = [0; 0];
guess(iphase).parameter    = [];

iphase = 2;
limits(iphase).time.min    = [t1 t2];
limits(iphase).time.max    = [t1 t2];
limits(iphase).state.min(1,:) = [rmin rmin rmin];
limits(iphase).state.max(1,:) = [rmax rmax rmax];
limits(iphase).state.min(2,:) = [rmin rmin rmin];
limits(iphase).state.max(2,:) = [rmax rmax rmax];
limits(iphase).state.min(3,:) = [rmin rmin rmin];
limits(iphase).state.max(3,:) = [rmax rmax rmax];
limits(iphase).state.min(4,:) = [vmin vmin vmin];
limits(iphase).state.max(4,:) = [vmax vmax vmax];
limits(iphase).state.min(5,:) = [vmin vmin vmin];
limits(iphase).state.max(5,:) = [vmax vmax vmax];
limits(iphase).state.min(6,:) = [vmin vmin vmin];
limits(iphase).state.max(6,:) = [vmax vmax vmax];
limits(iphase).state.min(7,:) = [m2f m2f m2f];
limits(iphase).state.max(7,:) = [m20  m20  m20];
limits(iphase).control.min(1,:) = -1;
limits(iphase).control.max(1,:) =  1;
limits(iphase).control.min(2,:) = -1;
limits(iphase).control.max(2,:) =  1;
limits(iphase).control.min(3,:) = -1;
limits(iphase).control.max(3,:) =  1;
limits(iphase).parameter.min    = [];
limits(iphase).parameter.max    = [];
limits(iphase).path.min    = 1;
limits(iphase).path.max    = 1;
guess(iphase).time = [t1; t2];
guess(iphase).state(:,1) = [r0(1); r0(1)];
guess(iphase).state(:,2) = [r0(2); r0(2)];
guess(iphase).state(:,3) = [r0(3); r0(3)];
guess(iphase).state(:,4) = [v0(1); v0(1)];
guess(iphase).state(:,5) = [v0(2); v0(2)];
guess(iphase).state(:,6) = [v0(3); v0(3)];
guess(iphase).state(:,7) = [m20; m2f];
guess(iphase).control(:,1) = [0; 0];
guess(iphase).control(:,2) = [1; 1];
guess(iphase).control(:,3) = [0; 0];
guess(iphase).parameter    = [];

iphase = 3;
limits(iphase).time.min = [t2 t3];
limits(iphase).time.max = [t2 t3];
limits(iphase).state.min(1,:) = [rmin rmin rmin];
limits(iphase).state.max(1,:) = [rmax rmax rmax];
limits(iphase).state.min(2,:) = [rmin rmin rmin];
limits(iphase).state.max(2,:) = [rmax rmax rmax];
limits(iphase).state.min(3,:) = [rmin rmin rmin];
limits(iphase).state.max(3,:) = [rmax rmax rmax];
limits(iphase).state.min(4,:) = [vmin vmin vmin];
limits(iphase).state.max(4,:) = [vmax vmax vmax];
limits(iphase).state.min(5,:) = [vmin vmin vmin];
limits(iphase).state.max(5,:) = [vmax vmax vmax];
limits(iphase).state.min(6,:) = [vmin vmin vmin];
limits(iphase).state.max(6,:) = [vmax vmax vmax];
limits(iphase).state.min(7,:) = [m3f m3f m3f];
limits(iphase).state.max(7,:) = [m30  m30  m30];
limits(iphase).control.min(1,:) = -1;
limits(iphase).control.max(1,:) =  1;
limits(iphase).control.min(2,:) = -1;
limits(iphase).control.max(2,:) =  1;
limits(iphase).control.min(3,:) = -1;
limits(iphase).control.max(3,:) =  1;
limits(iphase).parameter.min    = [];
limits(iphase).parameter.max    = [];
limits(iphase).path.min    = 1;
limits(iphase).path.max    = 1;
guess(iphase).time = [t2; t3];
guess(iphase).state(:,1) = [rout(1); rout(1)];
guess(iphase).state(:,2) = [rout(2); rout(2)];
guess(iphase).state(:,3) = [rout(3); rout(3)];
guess(iphase).state(:,4) = [vout(1); vout(1)];
guess(iphase).state(:,5) = [vout(2); vout(2)];
guess(iphase).state(:,6) = [vout(3); vout(3)];
guess(iphase).state(:,7) = [m30; m3f];
guess(iphase).control(:,1) = [0; 0];
guess(iphase).control(:,2) = [1; 1];
guess(iphase).control(:,3) = [0; 0];
guess(iphase).parameter    = [];

iphase = 4;
limits(iphase).time.min = [t3 t3];
limits(iphase).time.max = [t3 t4];
limits(iphase).state.min(1,:) = [rmin rmin rmin];
limits(iphase).state.max(1,:) = [rmax rmax rmax];
limits(iphase).state.min(2,:) = [rmin rmin rmin];
limits(iphase).state.max(2,:) = [rmax rmax rmax];
limits(iphase).state.min(3,:) = [rmin rmin rmin];
limits(iphase).state.max(3,:) = [rmax rmax rmax];
limits(iphase).state.min(4,:) = [vmin vmin vmin];
limits(iphase).state.max(4,:) = [vmax vmax vmax];
limits(iphase).state.min(5,:) = [vmin vmin vmin];
limits(iphase).state.max(5,:) = [vmax vmax vmax];
limits(iphase).state.min(6,:) = [vmin vmin vmin];
limits(iphase).state.max(6,:) = [vmax vmax vmax];
limits(iphase).state.min(7,:) = [m4f m4f m4f];
limits(iphase).state.max(7,:) = [m40 m40 m40];
limits(iphase).control.min(1,:) = -1;
limits(iphase).control.max(1,:) =  1;
limits(iphase).control.min(2,:) = -1;
limits(iphase).control.max(2,:) =  1;
limits(iphase).control.min(3,:) = -1;
limits(iphase).control.max(3,:) =  1;
limits(iphase).parameter.min = [];
limits(iphase).parameter.max = [];
limits(iphase).path.min      = 1;
limits(iphase).path.max      = 1;
limits(iphase).event.min     = [af; ef; incf; Omf; omf];
limits(iphase).event.max     = [af; ef; incf; Omf; omf];
guess(iphase).time    = [t3; t4];
guess(iphase).state(:,1) = [rout(1) rout(1)];
guess(iphase).state(:,2) = [rout(2) rout(2)];
guess(iphase).state(:,3) = [rout(3) rout(3)];
guess(iphase).state(:,4) = [vout(1) vout(1)];
guess(iphase).state(:,5) = [vout(2) vout(2)];
guess(iphase).state(:,6) = [vout(3) vout(3)];
guess(iphase).state(:,7) = [m40; m4f];
guess(iphase).control(:,1) = [0; 0];
guess(iphase).control(:,2) = [1; 1];
guess(iphase).control(:,3) = [0; 0];
guess(iphase).parameter    = [];

ipair = 1; % First pair of phases to link
linkages(ipair).left.phase = 1;
linkages(ipair).right.phase = 2;
linkages(ipair).min = [0; 0; 0; 0; 0; 0; -6*m_dry_srb];
linkages(ipair).max = [0; 0; 0; 0; 0; 0; -6*m_dry_srb];

ipair = 2; % Second pair of phases to link
linkages(ipair).left.phase = 2;
linkages(ipair).right.phase = 3;

linkages(ipair).min = [0; 0; 0; 0; 0; 0; -3*m_dry_srb];
linkages(ipair).max = [0; 0; 0; 0; 0; 0; -3*m_dry_srb];

ipair = 3; % Third pair of phases to link
linkages(ipair).left.phase = 3;
linkages(ipair).right.phase = 4;
linkages(ipair).min = [0; 0; 0; 0; 0; 0; -m_dry_first];
linkages(ipair).max = [0; 0; 0; 0; 0; 0; -m_dry_first];

setup.name = 'Launch-Vehicle-Ascent';
setup.funcs.cost = 'launchCost';
setup.funcs.dae = 'launchDae';
setup.funcs.event = 'launchEvent';
setup.funcs.link = 'launchLink';
setup.derivatives = 'finite-difference';
setup.checkDerivatives = 0;
setup.limits = limits;
setup.guess = guess;
setup.linkages = linkages;
setup.autoscale = 'off';

if isequal(setup.derivatives,'automatic-intlab'),
    CONSTANTS.derivatives = 'automatic-intlab';
else
    CONSTANTS.derivatives = [];
end;
setup.mesh.tolerance = 1e-6;
setup.mesh.iteration = 10;
setup.mesh.nodesPerInterval.min = 4;
setup.mesh.nodesPerInterval.max = 12;

tic
[output,gpopsHistory] = gpops(setup);
toc
solutionPlot= output.solutionPlot;
solution = output.solution;

