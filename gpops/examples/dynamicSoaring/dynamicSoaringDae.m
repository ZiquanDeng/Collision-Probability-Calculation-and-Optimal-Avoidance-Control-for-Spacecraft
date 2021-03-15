%-------------------------------------%
% BEGIN: function dynamicSoaringDae.m %
%-------------------------------------%
function [dae] = dynamicSoaringDae(soldae);

global CONSTANTS
t = soldae.time;
s = soldae.state;
u = soldae.control;
p = soldae.parameter;

x = s(:,1);
y = s(:,2);
z = s(:,3);
v = s(:,4);
r = s(:,5);
psi = s(:,6);
CL=u(:,1);
phi=u(:,2);
beta=p(1);

sinr = sin(r);
cosr = cos(r);
sinpsi = sin(psi);
cospsi = cos(psi);
sinphi = sin(phi);
cosphi = cos(phi);

rho=CONSTANTS.rho;
S=CONSTANTS.S;
CD0=CONSTANTS.CD0;
K=CONSTANTS.K;
g=CONSTANTS.g;
m=CONSTANTS.m;
ts=CONSTANTS.ts;
vs=CONSTANTS.vs;
ls=CONSTANTS.ls;
wx=(beta*ls*z+CONSTANTS.W0)/vs;
DWxDt=beta*vs*v.*sinr;


vcosr = v.*cosr;
DWxDtsinpsi = DWxDt.*sinpsi;

xdot=vcosr.*sinpsi+wx;
ydot=vcosr.*cospsi;
zdot=v.*sinr;
term1 = (rho*S*vs*ts/2/m);
term2 = ts/vs;
term3 = g*term2;
CLsq = CL.^2;
vsq = v.^2;
vdot=-term1*(CD0+K*CLsq).*vsq-(term3)*sinr-term2*DWxDtsinpsi.*cosr;
rdot= term1*CL.*v.*cosphi-(term3)*cosr./v+term2*DWxDtsinpsi.*sinr./v;
psidot=(term1*CL.*v.*sinphi-term2*DWxDt.*cospsi./v)./cosr;

% ng = load factor = L/(mg)
ngconstant = (0.5*CONSTANTS.rho*CONSTANTS.S/CONSTANTS.m/CONSTANTS.g);
ng=ngconstant.*CL.*(vs*v).^2;

dae = [xdot ydot zdot vdot rdot psidot ng];

%-----------------------------------%
% END: function dynamicSoaringDae.m %
%-----------------------------------%
