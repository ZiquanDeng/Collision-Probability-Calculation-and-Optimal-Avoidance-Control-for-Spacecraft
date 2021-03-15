function lagrange = gpopsSolutionInterpolation(setup,iphase)

nodesPerIntervalcum = cumsum(setup.limits(iphase).nodesPerInterval)+1;
indices = [1 nodesPerIntervalcum];
tau = [setup.ps(iphase).Points;1];
for seg = 1:size(setup.limits(iphase).nodesPerInterval,2)
  if (seg == 1)
    lagrange.T = [];
    lagrange.S = [];
    lagrange.C = [];
  end
  n = setup.limits(iphase).nodesPerInterval(seg);
  nstates = setup.sizes(iphase,1);
  ncontrols = setup.sizes(iphase,2);
  istart = indices(seg);
  ifinish = indices(seg+1);
  timestart  = tau(istart);
  timefinish = tau(ifinish);
  timecurrent = tau(istart:ifinish);
  statecurrent = setup.solution(iphase).state(istart:ifinish,:);
  if ~isempty(setup.solution(iphase).control),
    controlcurrent = setup.solution(iphase).control(istart:ifinish-1,:);
  end;
  costatecurrent = setup.solution(iphase).costate(istart:ifinish-1,:);
  Ttemp = [(((gpopsGetLGR(n+1)+1))*(timefinish-timestart)/2)+timestart;setup.limits(iphase).meshPoints(seg+1)];
  Stemp = zeros(size(Ttemp,1),setup.sizes(iphase,1));
  for i=1:nstates
    Stemp(:,i) = gpopsBaryLag([timecurrent,statecurrent(:,i)],Ttemp);
  end
  if ~isempty(setup.solution(iphase).control)
    Ctemp = zeros(size(Ttemp,1)-1,setup.sizes(iphase,2));
    for i=1:ncontrols
      Ctemp(:,i) = gpopsBaryLag([timecurrent(1:end-1,:),controlcurrent(:,i)],Ttemp(1:end-1));
    end
  else
    Ctemp = [];
  end    
  lagrange.T = [lagrange.T;Ttemp(1:end-1,:)];
  lagrange.S = [lagrange.S;Stemp(1:end-1,:)];
  lagrange.C = [lagrange.C;Ctemp(1:end,:)];
end
lagrange.T = [lagrange.T;1];
lagrange.S = [lagrange.S;setup.solution(iphase).state(end,:)];
