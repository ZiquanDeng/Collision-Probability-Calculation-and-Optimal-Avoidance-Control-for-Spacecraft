function [setup,endalgorithm] = gpopsRefineMesh(setup,tol)
%-------------------------------------------------------------------%
% This function refines the mesh on which the current solution to   %
% the NLP was obtained.  The mesh for the next execution of the NLP %
% solver is generated by determining if the error tolerance is      %
% satisfied in each segment of the problem.  If the error tolerance %
% in a mesh interval is satisfied, the mesh interval is left        %
% unchanged. If the mesh tolerance in a mesh interval is not        %
% satisfied, then either the degree of the polynomial is increased  %
% or the mesh interval is divided into subintervals.  The           %
% choice to increase the degree of the polynomial or subdivide the  %
% mesh interval is determined by the function GPOPSMODIFYSEGMENT.   %
%-------------------------------------------------------------------%

solution = setup.solution;
nphases = size(setup.limits,2);
endalgorithm = zeros(length(setup.limits),1);
for iphase = 1:nphases
  fprintf('\n');
  fprintf('\n');
  switch iphase<10
   case true
    fprintf(' _________________________________________________________________________________________\n');
    fprintf('|                                                                                         |\n');
    fprintf('|>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ANALYSIS OF MESH IN PHASE %d <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|\n',iphase);
    fprintf('|_________________________________________________________________________________________|\n');
   otherwise
    fprintf(' _________________________________________________________________________________________\n');
    fprintf('|                                                                                         |\n');
    fprintf('|>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ANALYSIS OF MESH IN PHASE %d <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|\n',iphase);
    fprintf('|_________________________________________________________________________________________|\n');
  end
  fprintf('\n');
  fprintf('\n');
  
  %--------------------------------------------------------------------%
  % Compute the relative error of the solution in the current phase    %
  % obtained by solving the NLP.                                       %
  %--------------------------------------------------------------------%
  relativeError = gpopsSolutionError(setup,iphase);
  %--------------------------------------------------------------------%
  % Now that the error has been computed, each mesh interval in the    %
  % solution must be traversed to see if the error exceeds the error   %
  % tolerance.  If the error exceeds the error tolerance, then either  %
  % the degree of the approximating polynomial in the segment is       %
  % increased or the mesh interval is subdivided.                      %
  %--------------------------------------------------------------------%
  numseg = size(setup.limits(iphase).nodesPerInterval,2);
  nodesPerIntervaluse = setup.limits(iphase).nodesPerInterval+1;
  daeSatisfied = zeros(numseg,1);
  stateScales = setup.stateScales{iphase};
  %--------------------------------------------------------------------%
  % Create a structure of length NUMSEG that contains a temporary mesh %
  %--------------------------------------------------------------------%
  TempMesh(1).nodesPerInterval = [];
  TempMesh(1).meshPoints = [];
  TempMesh(2:numseg) = TempMesh(1);
  for seg=1:numseg
    indexStart = sum(nodesPerIntervaluse(1:seg-1))+1;
    indexEnd   = sum(nodesPerIntervaluse(1:seg))+1;
    numIndices = length(indexStart:indexEnd);
    intervalError = relativeError(indexStart:indexEnd,:);
    nstates = setup.sizes(iphase,1);
    ErrorLessThanTolerance = (intervalError<=tol);
    daeSatisfied = all(ErrorLessThanTolerance(:));
    switch daeSatisfied
     case true
      % if daeSatisfied,
      %-----------------------------------------------------%
      % Do not modify the current mesh interval because the %
      % differential equations are satisfied.               %
      %-----------------------------------------------------%
      fprintf('Error Tolerance Satisfied in Segment %d of Phase %d\n',seg, iphase); 
      if ~exist('exitflagTEMP','var')
        exitflagTEMP = 0;
      end
     case false
      % else
      %--------------------------------------------------------------%
      % Modify the current mesh interval because the error tolerance %
      % on the differential equations is not satisfied.              %   
      %--------------------------------------------------------------%
      fprintf('Error Tolerance Not Satisfied in Segment %d of Phase %d\n',seg, iphase);
      [TempMesh,exitflagTEMP] = gpopsModifySegment(setup,iphase,intervalError,tol,seg,TempMesh);
    end
  end
  meshPointsPrev = setup.limits(iphase).meshPoints;
  nodesPerIntervalPrev = setup.limits(iphase).nodesPerInterval;
  meshPoints = [-1];
  nodesPerInterval = [];
  switch isempty(TempMesh)
   case false
    for seg = 1:(length(meshPointsPrev)-1)
      switch isempty(TempMesh(seg).meshPoints)
       case false
        meshPoints = [meshPoints,TempMesh(seg).meshPoints];
       case true
        meshPoints = [meshPoints,meshPointsPrev(seg+1)];
      end
      switch isempty(TempMesh(seg).nodesPerInterval)
       case false
        nodesPerInterval = [nodesPerInterval,TempMesh(seg).nodesPerInterval];    
       case true
        nodesPerInterval = [nodesPerInterval,nodesPerIntervalPrev(seg)];     
      end
    end
   case true
    meshPoints = meshPointsPrev;
    nodesPerInterval = nodesPerIntervalPrev;
  end
  
  if isequal(meshPoints,meshPointsPrev) && isequal(nodesPerInterval,nodesPerIntervalPrev) && isequal(exitflagTEMP,0)
    endalgorithm(iphase) = 1;
  elseif isequal(meshPoints,meshPointsPrev) && isequal(nodesPerInterval,nodesPerIntervalPrev) && isequal(exitflagTEMP,2)
    endalgorithm(iphase) = 2;
  else
    endalgorithm(iphase) = 3;
  end
  
  setup.limits(iphase).meshPoints = meshPoints;
  setup.limits(iphase).nodesPerInterval = nodesPerInterval;
  if setup.printoff == 0; 
    fprintf('\n'); 
  end
  
end % for iphase=1:numphases
