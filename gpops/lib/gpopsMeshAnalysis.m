function [setup,endalgorithm] = gpopsMeshAnalysis(tol,setup)
%------------------------------------------------------------------%
% This function determines the difference between two solution     %
% iterates over each segment.  Then, it is determined if           %
% a) The segment is accurate                                       %
% b) The segment needs more nodes                                  %
% c) The segment needs to be further segmented                     %
% A new setup structure is then created to resolve the problem     %
% over a new mesh                                                  %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

if 0,
solution = setup.solution;
nphases = size(setup.limits,2);
endalgorithm = zeros(length(setup.limits),1);
for iphase = 1:nphases
  fprintf('\n');
  % error = gpopsSolutionError(setup,iphase);
  % [setup,endalgorithm(iphase)] = gpopsCreateNewMesh(setup,iphase,error,tol);
  [setup,endalgorithm(iphase)] = gpopsCreateNewMesh(setup,iphase,tol);
end
if setup.printoff == 0; 
  fprintf('\n'); 
end
else
  [setup,endalgorithm] = gpopsCreateNewMesh(setup,tol);
end;

