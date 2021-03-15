function [TempMesh,exitflag] = gpopsModifySegment(setup,iphase,errorcurr,meshTolerance,seg,TempMesh)
%-------------------------------------------------------------------%
% This function modifies a mesh interval on the current mesh.  The  %
% decision to modify the mesh is made based on whether the error    %
% tolerance in the mesh interval is satisfied.  If the error is     %
% less than the error tolerance, then the mesh interval left        %
% unchanged.  Otherwise, either the degree of the approximating     %
% polynomial is increased or the mesh interval is divided into      %
% subintervals.  The degree of the approximating polynomial is      %
% increased if the ratio of the maximum curvature to the mean       %
% curvature is below the value RATIO.  The mesh interval is divided %
% if the ratio of the maximum curvature to the mean curvature is    %
% above the value RATIO.                                            %
%-------------------------------------------------------------------%

nodeshift = setup.mesh.nodesPerInterval.min;
splitmult = setup.mesh.splitmult;
nodesMin = setup.mesh.nodesPerInterval.min;
nodesMax = setup.mesh.nodesPerInterval.max;
nodesPerInterval = setup.limits(iphase).nodesPerInterval;
nodesCurr = nodesPerInterval(seg);
meshPoints = setup.limits(iphase).meshPoints;
state = setup.solution(iphase).state;
control = setup.solution(iphase).control;
maxes = max(abs(errorcurr));
[maxerror,whichstate]=max(maxes);
s = [setup.ps(iphase).Points; 1];
%-------------------------------------------%
% Compute the log10 of the ratio of the     %
% maximum error and the accuracy tolerance. %
%-------------------------------------------%
logdiff = ceil(log10(maxerror/meshTolerance));
%------------------------------------------------%
% if the value of LOGDIFF is less than zero, set %
% LOGDIFF to zero to ensure that the number of   %
% LGR points in the segment is not decreases.    %
%------------------------------------------------%
if logdiff < 0
  logdiff = 0;
end
%-----------------------------------------------------%
% If the algorithm decides that the segment must be   %
% divided, NUMBREAKS is the number of new intervals   %
% that are to be created.  In this process, it is     %
% necessary that the number of segment breaks be at   %
% least TWO.  Thus, numbreaks is NEVER less than TWO. %
%-----------------------------------------------------%
numbreaks = max(2,ceil(splitmult*logdiff));
%-------------------------------------------%
% Increase the number of LGR points in the  %
% current segment by the value of LOGDIFF.  %
%-------------------------------------------%
globalnodes = logdiff+nodesCurr;
if globalnodes > nodesMax,
  globalnodes = nodesMax;
end
if globalnodes < nodesMin,
  globalnodes = nodesMin;
end
magnitude = 0;
ratio = 2;

TempMesh(seg).meshPoints = [];
s0 = meshPoints(seg);
sf = meshPoints(seg+1);
nodesPerIntervalcum = cumsum(nodesPerInterval)+1;
indices = [1 nodesPerIntervalcum];
istart =indices(seg);
ifinish = indices(seg+1);
tcurr = s(istart:ifinish);
statecurr = state(istart:ifinish,:);
statecurrmaxerr = statecurr(:,whichstate);
%----------------------------------------------------%
% Compute the curvature in the current mesh interval %
%----------------------------------------------------%
curvatureRatio = gpopsCurvature(tcurr,statecurrmaxerr);
%-----------------------------------------------------------------%
% The code below performs the following check and analysis. First %
% it checks to see if the ratio of the maximum to the average     %
% curvature in the segment is less that the allowable value       %
% 'RATIO'. If this ratio is less than 'RATIO', then the algorithm %
% changes the degree of degree of the polynomial in the segment.  %
% If this ratio is greater than or equal to RATIO, then the       %
% algorithm divides the segment.                                  %
%-----------------------------------------------------------------%
curvatureCriterion = (curvatureRatio<ratio);
degreeCriterion1   = (nodesCurr~=nodesMax);
degreeCriterion2   = (nodesCurr==nodesMax);
if curvatureCriterion && degreeCriterion1,
  %---------------------------------------------------------------------%  
  % If this section of code is encountered, the error tolerance in the  %
  % mesh interval is NOT satisfied, but the curvature ratio is less     %
  % than the threshold given in RATIO.  As a result, the degree of the  %
  % approximating polynomial in the mesh interval is increased.         %
  %---------------------------------------------------------------------%  
  if setup.printoff == 0
    fprintf('Curvature Ratio Satisfied in Segment %d for phase %d\n',seg,iphase);
    fprintf('Curvature Ratio: %d\n',curvatureRatio);
    fprintf('Increased Poly. Degree from %d to %d\n',nodesCurr,globalnodes);
    fprintf('Logdiff %d\n',logdiff);
  end
  TempMesh(seg).nodesPerInterval = globalnodes;
  TempMesh(seg).meshPoints = [sf];
elseif ((curvatureCriterion && degreeCriterion2) || ~curvatureCriterion),
  %---------------------------------------------------------------------%  
  % If this section of code is encountered, neither the error tolerance %
  % in the mesh interval is satisfied NOR is the curvature ratio less   %
  % than the threshold given in RATIO.  As a result, the mesh interval  %
  % is divided into NUMBREAKS equally spaced mesh subintervals.         %
  %---------------------------------------------------------------------%  
  if setup.printoff == 0,
    fprintf('Curvature Ratio Not Satisfied in Segment %d for phase %d\n',seg, iphase);
  end
  limits = s0:(sf-s0)/numbreaks:sf;
  TempMesh(seg).nodesPerInterval = nodeshift*ones(1,numbreaks);
  TempMesh(seg).meshPoints = limits(2:end);
end
exitflag = 2;
