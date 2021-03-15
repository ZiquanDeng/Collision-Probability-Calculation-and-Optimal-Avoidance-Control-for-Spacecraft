function [setup,gpopsHistory] = gpopsSolveNLPandRefineMesh(setup)
%------------------------------------------------------------------%
% This function is the primary program that runs the main GPOPS    %
% driver and refines the mesh                                      %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David A. Benson              %
%------------------------------------------------------------------%
global igrid

%------------------------------------------------------------------%
% Check for any user-defined mesh parameters and set any undefined %
% parameters to their default values.  The DEFAULT values for the  %
% mesh refinement method are given as follows:                     %
%     setup.mesh.tolerance = 1e-3;                                 %
%     setup.mesh.iteration = 10;                                   %
%     setup.mesh.nodesPerInterval.min = 4;                         %
%     setup.mesh.nodesPerInterval.max = 12;                        %
%     setup.mesh.splitmult = 1.2;                                  %
%------------------------------------------------------------------%
switch isfield(setup,'mesh')
  case false
   setup.mesh.tolerance = 1e-3;
   setup.mesh.iteration = 10;
   setup.mesh.nodesPerInterval.min = 4;
   setup.mesh.nodesPerInterval.max = 12;
   setup.mesh.splitmult = 1.2;
 case true
  if ~isfield(setup.mesh,'tolerance') || isempty(setup.mesh.tolerance)
    setup.mesh.tolerance = 1e-3;
  end
  if ~isfield(setup.mesh,'iteration') || isempty(setup.mesh.iteration)
    setup.mesh.iteration = 10;
  end
  if ~isfield(setup.mesh,'splitmult') || isempty(setup.mesh.splitmult)
    setup.mesh.splitmult = 1.2;
  end
  switch isfield(setup.mesh,'nodesPerInterval')
    case false
     setup.mesh.nodesPerInterval.min = 4;
     setup.mesh.nodesPerInterval.max = 12;
   case true
    if ~isfield(setup.mesh.nodesPerInterval,'min') || isempty(setup.mesh.nodesPerInterval.min)
      setup.mesh.nodesPerInterval.min = 4;
    end
    if ~isfield(setup.mesh.nodesPerInterval,'max') || isempty(setup.mesh.nodesPerInterval.max)
      setup.mesh.nodesPerInterval.max = 12;
    end
  end
end

%------------------------------------------------------%
% Set the default grid for GPOPS if Not User-Specified %
%------------------------------------------------------%
defaultgrid = zeros(length(setup.limits),1);
for iphase = 1:length(setup.limits);
  if ~isfield(setup.limits(iphase),'meshPoints') || isempty(setup.limits(iphase).meshPoints) 
    if isfield(setup.limits(iphase),'nodesPerInterval') && ~isempty(setup.limits(iphase).nodesPerInterval)
      setup.limits(iphase).meshPoints = linspace(-1,1,length(setup.limits(iphase).nodesPerInterval)+1);
    else
      setup.limits(iphase).meshPoints = [-1 1];
    end        
    defaultgrid(iphase) = 1;
  elseif setup.limits(iphase).meshPoints(1) ~= -1 || setup.limits(iphase).meshPoints(end) ~= 1
    error('setup.limits(iphase).meshPoints must span -1 to +1 in phase %i',iphase)
  end
  if ~isfield(setup.limits(iphase),'nodesPerInterval') || isempty(setup.limits(iphase).nodesPerInterval)
    setup.limits(iphase).nodesPerInterval = repmat(20,1,length(setup.limits(iphase).meshPoints)-1);
    defaultgrid(iphase) = 1;
  end
  if length(setup.limits(iphase).meshPoints) ~= length(setup.limits(iphase).nodesPerInterval)+1
    error(['Number of nodesPerInterval must match number of mesh intervals in phase %i\n',...
           'i.e. length(setup.limits(iphase).meshPoints) == length(setup.limits(iphase).nodesPerInterval)+1'],iphase)
  end
end

if all(defaultgrid)
  fprintf('\n');
  fprintf('\n');
  fprintf(' _______________________________________________________________________________________________________\n');
  fprintf('|                                                                                                       |\n');
  fprintf('|>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INITIAL RUN WITH DEFAULT MESH <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|\n');
  fprintf('|_______________________________________________________________________________________________________|\n');
  fprintf('\n');
  fprintf('\n');
elseif any(defaultgrid)
  fprintf('\n');
  fprintf('\n');
  fprintf(' _______________________________________________________________________________________________________\n');
  fprintf('|                                                                                                       |\n');
  fprintf('|>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INITIAL RUN WITH USER-SPECIFIED MESH <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|\n');
  fprintf('|_______________________________________________________________________________________________________|\n');
  fprintf('\n');
  fprintf('\n');
else
  fprintf('\n');
  fprintf('\n');
  fprintf(' _______________________________________________________________________________________________________\n');
  fprintf('|                                                                                                       |\n');
  fprintf('|>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> INITIAL RUN WITH USER-SPECIFIED MESH <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|\n');
  fprintf('|_______________________________________________________________________________________________________|\n');
  fprintf('\n');
  fprintf('\n');
end

igrid = 1;
setup = gpopsSolveNLP(setup);
if nargout == 2
  gpopsHistory(igrid).output = setup;
end
tol = setup.mesh.tolerance;
iter = setup.mesh.iteration;
endalgorithm = 3*(ones(length(setup.limits),1));
while any(endalgorithm == 3) && (igrid<=iter)
  % [setup,endalgorithm] = gpopsMeshAnalysis(tol,setup);
  % [setup,endalgorithm] = gpopsCreateNewMesh(setup,tol);
  [setup,endalgorithm] = gpopsRefineMesh(setup,tol);
  if any(endalgorithm == 3)
    igrid = igrid + 1;
    fprintf('\n');
    fprintf('\n');
    if igrid<11,
      fprintf(' ___________________________________________________________________________________________\n');
      fprintf('|                                                                                           |\n');
      fprintf('|>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MESH REFINEMENT ITERATION %d <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|\n',igrid-1);
      fprintf('|___________________________________________________________________________________________|\n');
    else
      fprintf(' ____________________________________________________________________________________________\n');
      fprintf('|                                                                                            |\n');
      fprintf('|>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> MESH REFINEMENT ITERATION %d <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|\n',igrid-1);
      fprintf('|____________________________________________________________________________________________|\n');
    end
    fprintf('\n');
    fprintf('\n');
    setup.guess = setup.solution;
    [setup] = gpopsSolveNLP(setup);
    if nargout == 2
      gpopsHistory(igrid).output = setup;
    end
  end
end

if any(endalgorithm == 2)
  disp('                                                                                ');
  disp(' ________________________________________________________________________________|');
  disp('|                                                                                |');
  disp('| MESH REFINEMENT COMPLETED, BUT REQUESTED MESH TOLERANCE COULD NOT BE SATISFIED |');
  disp('|________________________________________________________________________________|');
  disp('                                                                                ');
end

if all(endalgorithm == 1)
  disp('                                                       ');
  disp(' _______________________________________________________');
  disp('|                                                       |');
  disp('| MESH REFINEMENT COMPLETED TO REQUESTED MESH TOLERANCE |');
  disp('|_______________________________________________________|');
  disp('                                                         ');
end

if igrid > iter
  disp('                                          ');
  disp(' _________________________________________');
  disp('|                                         |');
  disp('| MESH REFINEMENT ITERATION LIMIT REACHED |');
  disp('|_________________________________________|');
  disp('                                         ');
end

if igrid<10,
  disp(' ______________________________');
  disp('|                              |');
  disp(['| NUMBER OF SOLUTION MESHES: ',num2str(igrid),' |']);
  disp('|______________________________|');
else
  disp(' _______________________________');
  disp('|                               |');
  disp(['| NUMBER OF SOLUTION MESHES: ',num2str(igrid),' |']);
  disp('|_______________________________|');
end;  
disp('                                         ');
if igrid<10 
  disp(' _______________________________');
  disp('|                               |');
  disp(['| NUMBER OF MESH REFINEMENTS: ',num2str(igrid-1),' |']);
  disp('|_______________________________|');
else
  disp(' ________________________________');
  disp('|                                |');
  disp(['| NUMBER OF MESH REFINEMENTS: ',num2str(igrid-1),' |']);
  disp('|________________________________|');
end
disp('                                         ');
disp(' ___________________________________________________________________________________');
disp('|                                                                                   |');
disp('| SOLUTION AT DISCRETIZATION POINTS STORED IN "output.solution" OF OUTPUT STRUCTURE |');
disp('|___________________________________________________________________________________|');
disp('|                                                                                   |');
disp('|   solution.time      --> Array of Structures with Time in Each Phase              |');
disp('|   solution.state     --> Array of Structures with State in Each Phase             |');
disp('|   solution.control   --> Array of Structures with Control in Each Phase           |');
disp('|   solution.costate   --> Array of Structures with Costate in Each Phase           |');
disp('|   solution.parameter --> Array of Structures with Parameters in Each Phase        |');
disp('|___________________________________________________________________________________|');
disp('                                                                                    ');
disp(' ___________________________________________________________________________________');
disp('|                                                                                   |');
disp('|                SOLUTION FOR PLOTTING STORED in "output.solutionPlot"              |') 
disp('|___________________________________________________________________________________|');
disp(' ___________________________________________________________________________________');
disp('|                                                                                   |');
disp('|   solutionPlot.time      --> Array of Structures with Time in Each Phase          |');
disp('|   solutionPlot.state     --> Array of Structures with State in Each Phase         |');
disp('|   solutionPlot.control   --> Array of Structures with Control in Each Phase       |');
disp('|   solutionPlot.costate   --> Array of Structures with Costate in Each Phase       |');
disp('|   solutionPlot.parameter --> Array of Structures with Parameters in Each Phase    |');
disp('|___________________________________________________________________________________|');
disp('                                                                                    ');
