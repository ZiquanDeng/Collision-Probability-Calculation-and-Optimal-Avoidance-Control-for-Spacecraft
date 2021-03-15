function [Mayer,Lagrange,DMayer,DLagrange]=dynamicSoaringCost(solcost);

t0 = solcost.initial.time;
s0 = solcost.initial.state;
tf = solcost.terminal.time;
sf = solcost.terminal.state;
t  = solcost.time;
u = solcost.control;
s  = solcost.state;
p  = solcost.parameter;

% x = s(:,1);
% z = s(:,2);

Mayer = 7*p(1);
Lagrange = zeros(size(t));;

if nargout == 4
    DMayerdx0 = zeros(1,length(s0));
    DMayerdt0 = 0;
    DMayerdxf = zeros(1,length(sf));
    DMayerdtf = 0;
    DMayerdp  = 7;
    DMayer = [DMayerdx0 DMayerdt0 DMayerdxf DMayerdtf DMayerdp];
    DLagrangedx = zeros(length(t),length(s));
    DLagrangedu = zeros(length(t),length(u));
    DLagrangedt = zeros(length(t),1);
    DLagrangedp = zeros(length(t),length(p));
    DLagrange = [DLagrangedx DLagrangedu DLagrangedt DLagrangedp];
%    DLagrange =zeros(length(t),length(s)+length(u)+1+length(p));
 
end
