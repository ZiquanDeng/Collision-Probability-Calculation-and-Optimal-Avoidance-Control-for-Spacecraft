function [Sjac_I, Sjac_J, Sjac_V, Sconstant_I, Sconstant_J, Sconstant_V, Snonconstant_I, Snonconstant_J, Snonconstant_V] = ...
  gpopsPhaseSparsity(iphase,setup,phase_info,dependencies)
%------------------------------------------------------------------%
% Compute sparsity pattern for a single phase                      %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

nodes = phase_info(1);
sumnodes = nodes;
nstates = phase_info(2);
ncontrols = phase_info(3);
nparameters = phase_info(4);
npaths = phase_info(5);
nevents = phase_info(6);
ndiffeqs = nstates;
disc_pts = sumnodes+1;
numcons = ndiffeqs*(disc_pts-1)+npaths*sumnodes+nevents;
% find number of non-zero elements in dependencies (always including diagonal)
for i = 1:ndiffeqs
  dependencies(i,i) = 1; % set diag elements to one
end
nDependancies = nnz(dependencies);
% allocate memory for Jacobian
Sjac_I = zeros(nDependancies*sumnodes + 2*ndiffeqs*sumnodes + ...
  (ndiffeqs+npaths)*nparameters*sumnodes + nevents*(2*nstates+nparameters), 1);
Sjac_J = zeros(nDependancies*sumnodes + 2*ndiffeqs*sumnodes + ...
  (ndiffeqs+npaths)*nparameters*sumnodes + nevents*(2*nstates+nparameters), 1);
Sjac_V = zeros(nDependancies*sumnodes + 2*ndiffeqs*sumnodes + ...
  (ndiffeqs+npaths)*nparameters*sumnodes + nevents*(2*nstates+nparameters), 1);
% find non-zeros in off diagonal D matrix
[DiffMatOffDiag_I, DiffMatOffDiag_J, DiffMatOffDiag_V] = find(setup.ps(iphase).D-setup.ps(iphase).Ddiag);
nonZerosDiffMat = length(DiffMatOffDiag_I);
% find non-zeros in diagonal part ot D matrix
[DiffMatDiag_I,DiffMatDiag_J,DiffMatDiag_V] = find(setup.ps(iphase).Ddiag);
nonZerosDiffMatNonConstant = length(DiffMatDiag_I);
% allocate memory for constant derivatives
Sconstant_I = zeros(ndiffeqs*nonZerosDiffMat,1);
Sconstant_J = zeros(ndiffeqs*nonZerosDiffMat,1);
Sconstant_V = zeros(ndiffeqs*nonZerosDiffMat,1);
Snonconstant_I = zeros(ndiffeqs*nonZerosDiffMatNonConstant,1);
Snonconstant_J = zeros(ndiffeqs*nonZerosDiffMatNonConstant,1);
Snonconstant_V = zeros(ndiffeqs*nonZerosDiffMatNonConstant,1);
Sjac_rowShift = 0;
%
% differential eqns
%
for i=1:ndiffeqs;
    rowstart = (i-1)*(sumnodes);
    for j=1:nstates;
        colstart = (j-1)*disc_pts;
        if isequal(i,j),    
            % DdiagQuadBlock
            Sjac_rows = Sjac_rowShift + (1:sumnodes);
            Sjac_I(Sjac_rows) = rowstart + (1:sumnodes).';
            Sjac_J(Sjac_rows) = colstart + (1:sumnodes).';
            Sjac_V(Sjac_rows) = 1;
            Sjac_rowShift = Sjac_rowShift + sumnodes;
            % DiffMatOffDiag
            Sconstant_rows = (i-1)*nonZerosDiffMat + (1:nonZerosDiffMat);
            Sconstant_I(Sconstant_rows) = rowstart+DiffMatOffDiag_I;
            Sconstant_J(Sconstant_rows) = colstart+DiffMatOffDiag_J;
            Sconstant_V(Sconstant_rows) = DiffMatOffDiag_V;
            % DiffMatOffDiag
            Snonconstant_rows = (i-1)*nonZerosDiffMatNonConstant + (1:nonZerosDiffMatNonConstant);
            Snonconstant_I(Snonconstant_rows) = rowstart+DiffMatDiag_I;
            Snonconstant_J(Snonconstant_rows) = colstart+DiffMatDiag_J;
            Snonconstant_V(Snonconstant_rows) = DiffMatDiag_V;        
        else
            if isequal(dependencies(i,j),1), 
                % OffDiag_State
                Sjac_rows = Sjac_rowShift + (1:sumnodes);
                Sjac_I(Sjac_rows) = rowstart + (1:sumnodes).';
                Sjac_J(Sjac_rows) = colstart + (1:sumnodes).';
                Sjac_V(Sjac_rows) = 1;
                Sjac_rowShift = Sjac_rowShift + sumnodes;
            end;
        end;
    end;
    colshift = nstates*disc_pts;
    for j=1:ncontrols;
        colstart = colshift+(j-1)*sumnodes;       
        if isequal(dependencies(i,j+nstates),1)   
            % Diffeq_Block_Control
            Sjac_rows = Sjac_rowShift + (1:sumnodes);
            Sjac_I(Sjac_rows) = rowstart + (1:sumnodes).';
            Sjac_J(Sjac_rows) = colstart + (1:sumnodes).';
            Sjac_V(Sjac_rows) = 1;
            Sjac_rowShift = Sjac_rowShift + sumnodes;
        end;
    end;
end;
%
% path constraints
%
rowshift = ndiffeqs*(sumnodes);
for i=1:npaths;
    rowstart = rowshift+(i-1)*sumnodes;
    for j=1:nstates;
        colstart = (j-1)*disc_pts;
        if isequal(dependencies(i+ndiffeqs,j),1)  
            % Path_Block_State
            Sjac_rows = Sjac_rowShift + (1:sumnodes);
            Sjac_I(Sjac_rows) = rowstart + (1:sumnodes).';
            Sjac_J(Sjac_rows) = colstart + (1:sumnodes).';
            Sjac_V(Sjac_rows) = 1;
            Sjac_rowShift = Sjac_rowShift + sumnodes;
        end
    end;
    colshift = nstates*disc_pts;
    for j=1:ncontrols;
        colstart = colshift+(j-1)*sumnodes;
        if isequal(dependencies(i+ndiffeqs,j+nstates),1)  
            % Path_Block_Control
            Sjac_rows = Sjac_rowShift + (1:sumnodes);
            Sjac_I(Sjac_rows) = rowstart + (1:sumnodes).';
            Sjac_J(Sjac_rows) = colstart + (1:sumnodes).';
            Sjac_V(Sjac_rows) = 1;
            Sjac_rowShift = Sjac_rowShift + sumnodes;
        end;
    end;
end;
%
% event constraints
%
rowshift = ndiffeqs*(sumnodes)+npaths*sumnodes;
event_indices_init = 1:sumnodes+1:nstates*(sumnodes+1);
event_indices_term = sumnodes+1:sumnodes+1:nstates*(sumnodes+1);
for i = 1:nevents
    Sjac_rows = Sjac_rowShift + (1:2*nstates);
    Sjac_I(Sjac_rows) = rowshift+i;
    Sjac_J(Sjac_rows) = [event_indices_init.'; event_indices_term.'];
    Sjac_V(Sjac_rows) = 1;
    Sjac_rowShift = Sjac_rowShift + 2*nstates;
end
%
% to, tf
%
colshift = ndiffeqs*disc_pts+ncontrols*sumnodes;
Sjac_rows = Sjac_rowShift + (1:numcons);
Sjac_I(Sjac_rows) = (1:numcons).';
Sjac_J(Sjac_rows) = colshift+1;
Sjac_V(Sjac_rows) = 1;
Sjac_rowShift = Sjac_rowShift + numcons;
Sjac_rows = Sjac_rowShift + (1:numcons);
Sjac_I(Sjac_rows) = (1:numcons).';
Sjac_J(Sjac_rows) = colshift+2;
Sjac_V(Sjac_rows) = 1;
Sjac_rowShift = Sjac_rowShift + numcons;
%
% parameters
%
colshift = colshift+2;
for i = 1:nparameters
  Sjac_rows = Sjac_rowShift + (1:numcons);
  Sjac_I(Sjac_rows) = (1:numcons).';
  Sjac_J(Sjac_rows) = colshift+i;
  Sjac_V(Sjac_rows) = 1;
  Sjac_rowShift = Sjac_rowShift + numcons;
end

% numvars = nstates*disc_pts+ncontrols*sumnodes+nparameters+2;
% Sjac = sparse(Sjac_I,Sjac_J,Sjac_V,numcons,numvars);
% Sconstant = sparse(Sconstant_I,Sconstant_J,Sconstant_V,numcons,numvars);
