function [event Devent] = dynamicSoaringEvent(sol);

t0 = sol.initial.time;
x0 = sol.initial.state;
tf = sol.terminal.time;
xf = sol.terminal.state;
p  = sol.parameter;
%terminal constraints
%Dt=tf-t0;
event=[(xf(4)-x0(4));
    (xf(5)-x0(5));
    (xf(6)+2*pi-x0(6))];

if nargout == 2
    nevents = length(event);
    nstates = length(x0);
    nparameters = length(p);
    DEventdx0 = zeros(nevents,nstates);
    DEventdx0(1,4) = -1;
    DEventdx0(2,5) = -1;
    DEventdx0(3,6) = -1;
    DEventdt0 = zeros(nevents,1);
    DEventdxf = zeros(nevents,nstates);
    DEventdxf(1,4) = 1;
    DEventdxf(2,5) = 1;
    DEventdxf(3,6) = 1;
    DEventdtf = zeros(nevents,1);
    DEventdp  = zeros(length(event),nparameters);
    Devent = [DEventdx0 DEventdt0 DEventdxf DEventdtf DEventdp];
end


