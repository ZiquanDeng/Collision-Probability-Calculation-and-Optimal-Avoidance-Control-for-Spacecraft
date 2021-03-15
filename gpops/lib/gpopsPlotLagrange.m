function solutionPlot = gpopsPlotLagrange(setup)

solution = setup.solution;
solutionPlot = repmat(struct('time',[],'state',[],'control',[],'parameter',[],...
  'costate',[],'pathmult',[]),size(solution));
for iphase = 1:length(setup.limits)
  nodesPerIntervalcum = cumsum(setup.limits(iphase).nodesPerInterval)+1;
  indices = [1 nodesPerIntervalcum];
  for seg = 1:size(setup.limits(iphase).nodesPerInterval,2)
    n = setup.limits(iphase).nodesPerInterval(seg);        
    tau = [setup.ps(iphase).Points;1];
    istart = indices(seg);
    ifinish = indices(seg+1);
    timestart  = tau(istart);
    timefinish = tau(ifinish);
    timecurrent = tau(istart:ifinish);
    statecurrent = solution(iphase).state(istart:ifinish,:);
    costatecurrent = solution(iphase).costate(istart:ifinish-1,:);
    Ttemp = linspace(timestart,timefinish,5*n)';
    Stemp = zeros(size(Ttemp,1),setup.sizes(iphase,1));
    for i=1:setup.sizes(iphase,1)
      % Stemp(:,i) = gpopsLagrangeInterpolation(Ttemp,timecurrent,statecurrent(:,i));
      Stemp(:,i) = gpopsBaryLag([timecurrent,statecurrent(:,i)],Ttemp);
      % Stemp(:,i) = interp1(timecurrent,statecurrent(:,i),Ttemp,'pchip');
    end
    if ~isempty(solution(iphase).control)
      controlcurrent = solution(iphase).control(istart:ifinish-1,:);
      Ctemp = zeros(size(Ttemp,1),setup.sizes(iphase,2));
      for i=1:setup.sizes(iphase,2)
        % Ctemp(:,i) = gpopsLagrangeInterpolation(Ttemp,timecurrent(1:end-1,:),controlcurrent(:,i));
        % Ctemp(:,i) = gpopsBaryLag([timecurrent(1:end-1,:),controlcurrent(:,i)],Ttemp);
        Ctemp(:,i) = interp1(timecurrent(1:end-1,:),controlcurrent(:,i),Ttemp,'pchip');
      end
    else
      Ctemp = [];
    end
    CStemp = zeros(size(Ttemp,1),setup.sizes(iphase,1));
    for i=1:setup.sizes(iphase,1)
       % CStemp(:,i) = gpopsLagrangeInterpolation(Ttemp,timecurrent(1:end-1,:),costatecurrent(:,i)); 
       CStemp(:,i) = gpopsBaryLag([timecurrent(1:end-1,:),costatecurrent(:,i)],Ttemp);
       % CStemp(:,i) = interp1(timecurrent(1:end-1,:),costatecurrent(:,i),Ttemp,'pchip');
    end
    if ~isempty(solution(iphase).pathmult)
      pathmulticurrent = solution(iphase).pathmult(istart:ifinish-1,:);
      PMtemp = zeros(size(Ttemp,1),setup.sizes(iphase,4));
      for i=1:setup.sizes(iphase,4)
        % PMtemp(:,i) = gpopsLagrangeInterpolation(Ttemp,timecurrent(1:end-1,:),pathmulticurrent(:,i)); 
        % PMtemp(:,i) = gpopsBaryLag([timecurrent(1:end-1,:),pathmulticurrent(:,i)],Ttemp);
        PMtemp(:,i) = interp1(timecurrent(1:end-1,:),pathmulticurrent(:,i),Ttemp,'pchip');
      end
    else
    PMtemp = [];
    end
    if (seg == 1)
      solutionPlot(iphase).time = [solutionPlot(iphase).time;Ttemp(1:end,:)];
      solutionPlot(iphase).state = [solutionPlot(iphase).state;Stemp(1:end,:)];
      solutionPlot(iphase).control = [solutionPlot(iphase).control;Ctemp(1:end,:)];
      solutionPlot(iphase).costate = [solutionPlot(iphase).costate;CStemp(1:end,:)];
      solutionPlot(iphase).pathmult = [solutionPlot(iphase).pathmult;PMtemp(1:end,:)];
    else
      solutionPlot(iphase).time = [solutionPlot(iphase).time;Ttemp(2:end,:)];
      solutionPlot(iphase).state = [solutionPlot(iphase).state;Stemp(2:end,:)];
      solutionPlot(iphase).control = [solutionPlot(iphase).control;Ctemp(2:end,:)];
      solutionPlot(iphase).costate = [solutionPlot(iphase).costate;CStemp(2:end,:)];
      solutionPlot(iphase).pathmult = [solutionPlot(iphase).pathmult;PMtemp(2:end,:)];
    end
  end
  t0=solution(iphase).time(1);
  tf=solution(iphase).time(end);
  solutionPlot(iphase).time = (((solutionPlot(iphase).time+1)*.5)*(tf-t0))+t0;
  solutionPlot(iphase).parameter = solution(iphase).parameter;
end
