function [output,gpopsHistory] = gpops(setup)
%------------------------------------------------------------------%
%     GPOPS:  General Pseudopsectral Optimal Control Software      %
%------------------------------------------------------------------%
%  GPOPS is a MATLAB(R) program for solving multiple-phase optimal %
%  control problems. GPOPS uses an hp-adaptive version of the      %
%  Radau pseudospectral method (RPM).  The Radau pseudospectral    %
%  method employs collocation at the Legendre-Gauss-Radau points.  %
%  GPOPS is based entirely on mathematical methods that have been  %
%  published in the open literature.  Further information  about   %
%  the pseudospectral theory used in GPOPS and the applications    %
%  of this theory to various problems of interest can be found at  %
%  http://vdol.mae.ufl.edu.                                        % 
%                                                                  %
%------------------------------------------------------------------%
%  Authors of the GPOPS Code:                                      %
%    David A. Benson:      Draper Laboratory, Cambridge, MA        %
%    Anil V. Rao:          University of Florida, Gainesville, FL  %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David A. Benson              %
%------------------------------------------------------------------%
% For Licensing Information, Please See File LICENSE               %
% ---------------------------------------------------------------- %
%                                                                  %
%                SEE THE GPOPS MANUAL FOR USAGE                    %
%                                                                  %
%------------------------------------------------------------------%

%---------------------------------%
% Set print flag if it is missing %
%---------------------------------%
if ~isfield(setup,'printoff'); 
  setup.printoff = 0; 
end
%------------------------------------------------------%
% Print the GPOPS header to the MATLAB Command Window. %
%------------------------------------------------------%
gpopsSplash;
%-----------------------------------%
% Solve the NLP and refine the mesh %
%-----------------------------------%
if nargout == 2
  [output,gpopsHistory] = gpopsSolveNLPandRefineMesh(setup);
else
  output = gpopsSolveNLPandRefineMesh(setup);
end
%---------------%
% Add plot data %
%---------------%
output.solutionPlot = gpopsPlotLagrange(output);
%-----------------------------------------%
% Remove extra fields in output structure %
%-----------------------------------------%
output = gpopsClearFields(output);
%-------------------------------%
% Remove GPOPS global variables %
%-------------------------------%
clearvars -global igrid mysetup
