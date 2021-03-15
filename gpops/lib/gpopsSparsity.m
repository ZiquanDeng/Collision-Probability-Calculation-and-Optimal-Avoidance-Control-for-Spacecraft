function setup = gpopsSparsity(setup)
%------------------------------------------------------------------%
% Generate the sparsity pattern for a multiple-phase optimal       %
% control setup discretized using the Gauss pseudospectral         %
% method.                                                          %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

sizes = setup.sizes;
nodes = setup.nodes;
numphases = setup.numphases;
dependencies = setup.dependencies;
% find total non-zero elements in Jacobian and constant derivatives
nonZeros_Sjac = 0;
nonZeros_Sconstant = 0;
nonZeros_Snonconstant = 0;
for i=1:numphases;
  nstates = sizes(i,1);
  ncontrols = sizes(i,2);
  nparameters = sizes(i,3);
  npaths = sizes(i,4);
  nevents = sizes(i,5);
  % find number of non-zero elements in dependencies (always including diagonal)
  dependenciesTmp = dependencies{i};
  for j = 1:nstates
    dependenciesTmp(j,j) = 1; % set diag elements to one
  end
  nDependancies = nnz(dependenciesTmp);
  nonZeros_Sjac = nonZeros_Sjac + nDependancies*nodes(i) + 2*nstates*nodes(i) + ...
      (nstates+npaths)*nparameters*nodes(i) + nevents*(2*nstates+nparameters);
  nonZeros_Sconstant = nonZeros_Sconstant+nstates*nnz(setup.ps(i).D-setup.ps(i).Ddiag);
  nonZeros_Snonconstant = nonZeros_Snonconstant+nstates*nnz(setup.ps(i).Ddiag);
end
for ipair=1:setup.numlinkpairs;
  left_phase  = setup.linkages(ipair).left.phase;
  right_phase = setup.linkages(ipair).right.phase;
  nstates_left = sizes(left_phase,1);
  nparameters_left = sizes(left_phase,3);
  nstates_right = sizes(right_phase,1);
  nparameters_right = sizes(right_phase,3);
  numlinks = length(setup.linkages(ipair).min); 
  nonZeros_Sjac = nonZeros_Sjac + numlinks * (nstates_left + nparameters_left + nstates_right + nparameters_right);
end
% add objective row
nonZeros_Sjac = nonZeros_Sjac + setup.numvars; 
% allocate memory for Jacobian and constant derivatives
Sjac_I = zeros(nonZeros_Sjac,1);
Sjac_J = zeros(nonZeros_Sjac,1);
Sjac_V = zeros(nonZeros_Sjac,1);
Sjac_rowsShift = 0;
Sconstant_I = zeros(nonZeros_Sconstant,1);
Sconstant_J = zeros(nonZeros_Sconstant,1);
Sconstant_V = zeros(nonZeros_Sconstant,1);
Snonconstant_I = zeros(nonZeros_Snonconstant,1);
Snonconstant_J = zeros(nonZeros_Snonconstant,1);
Snonconstant_V = zeros(nonZeros_Snonconstant,1);
Sconstant_rowsShift = 0;
Snonconstant_rowsShift = 0;
rowshift = 0;
colshift = 0;
% ----------------------------------------
% Sparsity Pattern for each Phase
% ----------------------------------------
for i=1:numphases;
  nstates = sizes(i,1);
  ncontrols = sizes(i,2);
  nparameters = sizes(i,3);
  npaths = sizes(i,4);
  nevents = sizes(i,5);
  phase_info = [nodes(i); nstates; ncontrols; nparameters; npaths; nevents];
  [Sjac_phase_I, Sjac_phase_J, Sjac_phase_V, Sconstant_phase_I, Sconstant_phase_J, Sconstant_phase_V, Snonconstant_phase_I, Snonconstant_phase_J, Snonconstant_phase_V] ...
      = gpopsPhaseSparsity(i,setup,phase_info,dependencies{i});
  
  Sjac_rows = Sjac_rowsShift + (1:length(Sjac_phase_I));
  Sjac_I(Sjac_rows) = rowshift + 1 + Sjac_phase_I; % account for objective
  Sjac_J(Sjac_rows) = colshift + Sjac_phase_J;
  Sjac_V(Sjac_rows) = Sjac_phase_V;
  Sjac_rowsShift = Sjac_rowsShift + length(Sjac_phase_I);
  
  Sconstant_rows = Sconstant_rowsShift + (1:length(Sconstant_phase_I));
  Sconstant_I(Sconstant_rows) = rowshift + Sconstant_phase_I;
  Sconstant_J(Sconstant_rows) = colshift + Sconstant_phase_J;
  Sconstant_V(Sconstant_rows) = Sconstant_phase_V;    
  Sconstant_rowsShift = Sconstant_rowsShift + length(Sconstant_phase_I);
  
  Snonconstant_rows = Snonconstant_rowsShift + (1:length(Snonconstant_phase_I));
  Snonconstant_I(Snonconstant_rows) = rowshift + Snonconstant_phase_I;
  Snonconstant_J(Snonconstant_rows) = colshift + Snonconstant_phase_J;
  Snonconstant_V(Snonconstant_rows) = Snonconstant_phase_V;    
  Snonconstant_rowsShift = Snonconstant_rowsShift + length(Snonconstant_phase_I);
  
  numcons = nstates*nodes(i)+npaths*nodes(i)+nevents;
  numvars = nstates*(nodes(i)+1)+ncontrols*nodes(i)+nparameters+2;
  rowshift = rowshift+numcons;
  colshift = colshift+numvars;
end;
% ----------------------------------------
% Sparsity Pattern for Linkage Conditions
% ----------------------------------------
linkages = setup.linkages;
numlinkpairs = setup.numlinkpairs;
%numvars = setup.numvars;
variable_indices = setup.variable_indices;
SjacRows = setup.numnonlin - setup.numlinks + 1; % account for objective
for ipair=1:numlinkpairs;
  left_phase  = linkages(ipair).left.phase;
  right_phase = linkages(ipair).right.phase;
  nodes_left = nodes(left_phase);
  disc_left  = nodes_left+1;
  nstates_left = sizes(left_phase,1);
  ncontrols_left = sizes(left_phase,2);
  nparameters_left = sizes(left_phase,3);
  nodes_right = nodes(right_phase);
  disc_right  = nodes_right+1;
  nstates_right = sizes(right_phase,1);
  ncontrols_right = sizes(right_phase,2);
  nparameters_right = sizes(right_phase,3);
  nparameters_left_start = disc_left*nstates_left+nodes_left*ncontrols_left+3; %account for to,tf
  nparameters_left_finish = disc_left*nstates_left+nodes_left*ncontrols_left+2+nparameters_left; %account for to,tf
  parameter_indices_left = nparameters_left_start:nparameters_left_finish;
  indices_left = [disc_left:disc_left:nstates_left*disc_left, parameter_indices_left];
  nparameters_right_start  = disc_right*nstates_right+nodes_right*ncontrols_right+3; %account for to,tf
  nparameters_right_finish = disc_right*nstates_right+nodes_right*ncontrols_right+2+nparameters_right; %account for to,tf
  parameter_indices_right = nparameters_right_start:nparameters_right_finish;
  indices_right = [1:disc_right:nstates_right*disc_right, parameter_indices_right];
  indices_left_use = variable_indices{left_phase}(indices_left);
  indices_right_use = variable_indices{right_phase}(indices_right);
  numlinks = length(linkages(ipair).min);   
  for i = 1:numlinks
    Sjac_rows = Sjac_rowsShift + (1:length(indices_left_use));
    Sjac_I(Sjac_rows) = SjacRows+i; 
    Sjac_J(Sjac_rows) = indices_left_use.';
    Sjac_V(Sjac_rows) = 1;
    Sjac_rowsShift = Sjac_rowsShift + length(indices_left_use);
    Sjac_rows = Sjac_rowsShift + (1:length(indices_right_use));
    Sjac_I(Sjac_rows) = SjacRows+i; 
    Sjac_J(Sjac_rows) = indices_right_use.';
    Sjac_V(Sjac_rows) = 1;
    Sjac_rowsShift = Sjac_rowsShift + length(indices_right_use);
  end    
  SjacRows = SjacRows + numlinks;
end;
% ----------------------------------------
% Sparsity Pattern for Objective
% ----------------------------------------
Sjac_rows = Sjac_rowsShift + (1:setup.numvars);
Sjac_I(Sjac_rows) = 1; 
Sjac_J(Sjac_rows) = (1:setup.numvars).';
Sjac_V(Sjac_rows) = 1;

% store as sparse matrices
setup.sparsity_all = sparse(Sjac_I,Sjac_J,Sjac_V,1+setup.numnonlin+setup.numlin,setup.numvars);
setup.sparsity_constant = sparse(Sconstant_I,Sconstant_J,Sconstant_V,setup.numnonlin,setup.numvars);
setup.sparsity_nonconstant = sparse(Snonconstant_I,Snonconstant_J,Snonconstant_V,setup.numnonlin,setup.numvars);
[iRowsConstant,jColumnsConstant] = size(setup.sparsity_constant);
[iRowsnonConstant,jnonColumnsConstant] = size(setup.sparsity_nonconstant);
[iRowsAll,jColumnsAll] = size(setup.sparsity_all);
setup.sparsity = sparse(iRowsAll,jColumnsAll);
setup.sparsity(1,:) = setup.sparsity_all(1,:);
setup.sparsity(2:iRowsConstant+1,:) = sign(setup.sparsity_all(2:iRowsConstant+1,:)+setup.sparsity_constant);
