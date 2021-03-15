function setup = gpopsSolveNLP(setup)
%-------------------------------------------------------------------%
% This function solves the nonlinear programming problem (NLP) on a %
% mesh that has been provided.  If this is the first call to this   %
% function, either GPOPS provides a default mesh or the user has    %
% specified an initial mesh.  If is is not the first call to this   %
% function, the mesh is provided by the mesh refinement algorithm.  %
% Given a mesh, this function then sets up and solves the NLP using %
% any NLP solver for which GPOPS has been properly configured.      %
%-------------------------------------------------------------------%

%-------------------------------------------------------------------%
% The SNOPT MATLAB interface uses a fixed number of input arguments.%
% Consequently, user-defined inputs cannot be provided in the call  % 
% to SNOPT.  Because it is necessary to access information from the %
% SETUP structure in various functions, a copy of this structure is %
% made a global variable.  The copy of SETUP is called MYSETUP.     %
%-------------------------------------------------------------------%
global mysetup igrid
%------------------------------------------------------------------%
%        Clear any previous calls to the SNOPT mex file            %
%------------------------------------------------------------------%
clear snoptcmex;
%------------------------------------------------------------------%
%        Check for SNOPT installation                              %
%------------------------------------------------------------------%
if isempty(which('snoptcmex'))
  error('SNOPT not found or installed incorrectly')
end
%------------------------------------------------------------------%
%    Get the sizes in each phase of the optimal control problem    %
%------------------------------------------------------------------%
setup = gpopsGetSizes(setup);
%------------------------------------------------------------------%
%        Print a description of the optimal control problem        %
%------------------------------------------------------------------%
if igrid==1,
  gpopsPrint(setup);
end;
setup.igrid = igrid;
%------------------------------------------------------------------%
%    Get the bounds in each phase of the optimal control problem   %
%------------------------------------------------------------------%
setup = gpopsGetBounds(setup);
%------------------------------------------------------------------%
%               Get the initial guess for the NLP                  %
%------------------------------------------------------------------%
setup = gpopsGetGuess(setup);
%------------------------------------------------------------------%
%             Scale the nonlinear programming problem              %
%------------------------------------------------------------------%
setup = gpopsScaleNlp(setup);
colScales = setup.column_scales;
colShifts = setup.column_shifts;
rowScales = setup.row_scales;
varMin = setup.varbounds_min;
varMax = setup.varbounds_max;
conMin = setup.conbounds_min;
conMax = setup.conbounds_max;
%------------------------------------------------------------------%
%            Determine the sparsity pattern of the NLP             %
%------------------------------------------------------------------%
setup = gpopsSparsity(setup);
%------------------------------------------------------------------%
%           Lower and upper bounds on the NLP variables            %
%------------------------------------------------------------------%
xlow = varMin.*colScales+colShifts;
xupp = varMax.*colScales+colShifts;
%------------------------------------------------------------------%
% Setup diagonal matrices containing the scale factors computed by %
% the automatic scaling routine.                                   %
%------------------------------------------------------------------%
Dx = spdiags(colScales,0,length(colScales),length(colScales));
invDx = spdiags(1./colScales,0,length(colScales),length(colScales));
DF = spdiags(rowScales,0,length(rowScales),length(rowScales));
invDF = spdiags(1./rowScales,0,length(rowScales),length(rowScales));
setup.Dx=Dx; setup.invDx=invDx; setup.DF=DF; setup.invDF=invDF;
%------------------------------------------------------------------%
%            SCONSTANT contains the CONSTANT                       %
%            part of the linear derivattives                       %
%------------------------------------------------------------------%
Sconstant = setup.sparsity_constant;
Sconstant = DF*Sconstant*invDx;
%------------------------------------------------------------------%
%            SNONCONSTANT contains the NONCONSTANT                 %
%              part of the linear derivattives                     %
%------------------------------------------------------------------%
% Snonconstant = setup.sparsity_nonconstant;
%------------------------------------------------------------------%
%     Matrix containing coefficients of the linear constraints     %
%------------------------------------------------------------------%
Alinear = setup.Alinear;
Alinear = Alinear*invDx;
%------------------------------------------------------------------%
%  Modify the lower & upper bounds on the variables & constraints  %
%  depending upon whether or not the user has chosen to            %
%  automatically scale the NLP.                                    %
%------------------------------------------------------------------%
if isfield(setup,'autoscale') && isequal(setup.autoscale,'on'),
  disp('| Automatic Scaling Turned On                                                 |');
  %-------------------------------------------------- %
  % Give Warning for infinite bounds with autoscaling %
  %-------------------------------------------------- %
  if ~isempty(find(isinf(xlow)|isinf(xupp),1))
    disp('| WARNING!!! Automatic Scaling may not work with infinite bounds              |')
    disp(' ')
  end
else
  disp('| Automatic Scaling Turned Off                                                |');
end;
%------------------------------------------------------------------%
%        Lower and upper bounds on the linear constraints          %
%------------------------------------------------------------------%
Alinmin = setup.Alinmin;
Alinmax = setup.Alinmax;
%------------------------------------------------------------------%
%       Lower and upper bounds on the nonlinear constraints        %
%------------------------------------------------------------------%
clow = conMin;
cupp = conMax;
%------------------------------------------------------------------%
%                   Initial guess for the NLP                      %
%------------------------------------------------------------------%
xguess = setup.nlpGuess.*colScales+colShifts;
%------------------------------------------------------------------%
%             Find the row and column indices in S_ALL             %
%------------------------------------------------------------------%
[iGfun,jGvar]=find(setup.sparsity);
%------------------------------------------------------------------%
%             Sort the row and column indices by rows.             %
%------------------------------------------------------------------%
JJ = sortrows([iGfun jGvar],1); 
iGfun = JJ(:,1);
jGvar = JJ(:,2);
setup.iGfun = iGfun;
setup.jGvar = jGvar;
%------------------------------------------------------------------%
% ALINEAR_AUGMENTED is a matrix of the linear constraints plus the %
% constant derivatives.                                            %
%------------------------------------------------------------------%
Alinear_augmented = [zeros(1,setup.numvars); Sconstant; Alinear];
setup.Alinear_augmented = Alinear_augmented;
[iAfun,jAvar,AA] = find(Alinear_augmented);
%------------------------------------------------------------------%
%             Sort the row and column indices by rows              %
%------------------------------------------------------------------%
II = sortrows([iAfun jAvar AA],1);
iAfun = II(:,1);
jAvar = II(:,2);
AA    = II(:,3);
%------------------------------------------------------------------%
%        Free MATLAB memory by clearing unneeded variables         %
%------------------------------------------------------------------%
clear JJ II Alinear_augmented % Sconstant Alinear 

%------------------------------------------------------------------%
% This section sets the appropriate differentiation                %
% method for use with the NLP solver.                              %
%------------------------------------------------------------------%
% if isfield(setup,'derivatives'),
switch isfield(setup,'derivatives'),
  case true
   deropt = lower(setup.derivatives);
   if isequal(deropt,'automatic')
     %---------------------------------------------------%
     % Check whether Built-In automatic differentiation  %
     % is installed on the machine.                      %
     %---------------------------------------------------%
     if isempty(which('ad.m'))
       error('Built-in automatic differentiator not found or installed incorrectly')
     end
     %-----------------------------------------------------------%
     % If setup.derivatives equals 'automatic', then use the     %
     % built-in automatic differentiator for the computation of  %
     % derivatives.  When using BUILT-IN AUTOMATIC               %
     % DIFFERENTIATION, SNOPT will call the function             %
     % 'gpopsuserfunAD.                                          %
     %-----------------------------------------------------------%
     userfun = 'gpopsuserfunAD';
     snseti('Derivative Option',1);
     disp('| Objective Gradient Being Estimated via Built-In Automatic Differentiation   |');        
     disp('| Constraint Jacobian Being Estimated via Built-In Automatic Differentiation  |');
   elseif isequal(deropt,'automatic-intlab')
     %---------------------------------------------------%
     % Check whether INTLAB automatic differentiation    %
     % is installed on the machine.                      %
     %---------------------------------------------------%
     if isempty(which('gradientinit.m'))
       error('INTLAB not found or installed incorrectly')
     end
     %-----------------------------------------------------------%
     % If setup.derivatives equals EITHER 'automatic-INTLAB'     %
     % OR 'automatic', then use INTLAB for computing derivatives %
     % When using INTLAB AUTOMATIC DIFFERENTIATION, SNOPT will   %
     % function 'gpopsuserfunADINT.m'                            %
     %-----------------------------------------------------------%
     userfun = 'gpopsuserfunADINT';
     snseti('Derivative Option',1);
     disp('| Objective Gradient Being Estimated via INTLAB Automatic Differentiation     |');        
     disp('| Constraint Jacobian Being Estimated via INTLAB Automatic Differentiation    |');
   elseif isequal(deropt,'finite-difference'),
     %-----------------------------------------------------------%
     % When using NUMERICAL DIFFERENTIATION, SNOPT will call the %
     % function 'gpopsuserfunFD.m'                               %
     %-----------------------------------------------------------% 
     userfun = 'gpopsuserfunFD';
     snseti('Derivative Option',1);
     disp('| Objective Function Gradient Being Estimated via Sparse Finite-Differencing  |');        
     disp('| Constraint Jacobian Being Estimated via Sparse Finite-Differencing          |');
   elseif isequal(deropt,'complex'),
     %-----------------------------------------------------------%
     % When using COMPLEX-STEP DIFFERENTIATION, SNOPT will call  %
     % the function 'gpopsuserfunCS.m'                           %
     %-----------------------------------------------------------%
     userfun = 'gpopsuserfunCS';
     snseti('Derivative Option',1);
     setup.hpert = 1e-20;
     setup.deltaxmat = sqrt(-1)*setup.hpert*speye(setup.numvars);
     setup.Jaczeros = zeros(setup.numnonlin+setup.numlin+1,setup.numvars);
     disp('| Objective Gradient Being Estimated via Complex-Step Approximation           |');
     disp('| Constraint Jacobian Being Estimated via Complex-Step Approximation          |');
   elseif isequal(deropt,'analytic'),
     %-----------------------------------------------------------%
     % When using ANALYTIC DIFFERENTIATION, SNOPT will call the  %
     % function 'gpopsuserfunAN.m'                               %
     %-----------------------------------------------------------%
     userfun = 'gpopsuserfunAN';
     snseti('Derivative Option',1);
     disp('| Objective Gradient Being Estimated via Analytic Differentiation             |');        
     disp('| Constraint Jacobian Being Estimated via Analytic Differentiation            |');
     if isfield(setup,'checkDerivatives') && setup.checkDerivatives == 1
       %---------------------------------%
       % Check user supplied derivatives %
       %---------------------------------%
       gpopsCheckDerivatives(setup);
     end
   else
     error(['Unknown derivative option "%s" in setup.derivatives\n',...
            'Valid options are:',...
            '\n\tautomatic         \t-> GPOPS built-in Automatic Differentiation',...
            '\n\tautomatic-INTLAB  \t-> INTLAB Automatic Differentiation',...
            '\n\tfinite-difference \t-> Sparse Finite-Differencing',...
            '\n\tcomplex           \t-> Complex-Step Differentiation',...
            '\n\tanalytic          \t-> User Defined Analytic Differentiation\n'],...
           setup.derivatives);
   end;
% else
 case false
  %-----------------------------------------------------------%
  % When using FINITE-DIFFERENCE APPROXIMATIONS, SNOPT will   %
  % call the function 'gpopsuserfunFD.m'                      %
  %-----------------------------------------------------------%
  userfun = 'gpopsuserfunFD';
  % Set SNOPT Derivative Option==1 because finite-differencing
  % is being applied sparsely to optimal control problem.  Thus,
  % SNOPT thinks the derivative is analytic even though it is only
  % an approximation.   
  snseti('Derivative Option',1);
  disp('Objective Gradient Being Computed via Finite Differencing');
  disp('Constraint Jacobian Being Computed via Finite Differencing');
end;
fprintf('|_____________________________________________________________________________|\n');
fprintf('                                                                              \n');
%---------------------------------------------%
% NUMLIN    = Number of linear constraints    %
% NUMNONLIN = Number of nonlinear constraints %
%---------------------------------------------%
numlin    = setup.numlin;
numnonlin = setup.numnonlin;
%---------------------------------------------------------------%
% Pre-allocate a bunch of vectors for use in the user function. %
% This pre-allocation is done to reduce execution time.         %
%---------------------------------------------------------------%
setup.initlincons = zeros(numlin,size(xguess,2));
%----------------------------------------------------------%
% Set up settings that are used by SNOPT for every problem %
%----------------------------------------------------------%
snprint('snoptmain.out');         % Name of SNOPT Print File
snseti('Timing level',3);         % Print Execution Time to File
snset('Hessian Limited Memory');  % Choose Hessian Type
snset('LU Partial Pivoting');     % Choose Pivoting Method
snseti('Verify Level',-1);        % Derivative Verification Level
snset('Cold Start');
if isfield(setup,'maxIterations'),
  if isa(setup.maxIterations,'double'),
    maxIters = setup.maxIterations;
  else
    error('field <maxIterations> must be a double');
  end;
else
  maxIters = 100000;
end;
snseti('Iteration Limit',10*maxIters); % Iteration Limit
snseti('Major Iterations Limit',maxIters); % Major Iteration Limit
snseti('Minor Iterations Limit',maxIters); % Minor Iteration Limit
if isfield(setup,'tolerances');
  if isequal(length(setup.tolerances),2),
    if ~isempty(setup.tolerances(1)),
      snsetr('Major Optimality Tolerance',setup.tolerances(1));
      snsetr('Major Feasibility Tolerance',setup.tolerances(2));
    else
      snsetr('Major Optimality Tolerance',setup.tolerances(2));
    end;
  elseif isequal(length(setup.tolerances),1),
    snsetr('Major Optimality Tolerance',setup.tolerances(1));
  end;
else
  if isfield(setup,'autoscale'),
    if isequal(setup.autoscale,'on'),
      setup.tolerances(1) = 1e-6;
      setup.tolerances(2) = 2e-6;
      snsetr('Major Optimality Tolerance',setup.tolerances(1));
      snsetr('Major Feasibility Tolerance',setup.tolerances(2));
    else
      setup.tolerances(1) = 1e-6;
      setup.tolerances(2) = 2e-6;
      snsetr('Major Optimality Tolerance',setup.tolerances(1));
      snsetr('Major Feasibility Tolerance',setup.tolerances(2));
    end;
  else
    setup.autoscale='on';
    setup.tolerances(1) = 1e-6;
    setup.tolerances(2) = 2e-6;
    snsetr('Major Optimality Tolerance',setup.tolerances(1));
    snsetr('Major Feasibility Tolerance',setup.tolerances(2));
  end;
end;
%------------------------------------------------------------------%
% Set the lower and upper limits on the constraints and objective  %
% function.  The following assumptions are made:                   %
%  Row 1 is the objective row                                      %
%  Rows 2 through numnonlin+1 are the nonlinear constraints        %
%  Rows numnonlin+2 to the end are the linear constraints          %
%------------------------------------------------------------------%
nonlinearBoundsShift = DF*Sconstant*invDx*colShifts;
clow = DF*clow+nonlinearBoundsShift;
cupp = DF*cupp+nonlinearBoundsShift;
linearBoundsShift = Alinear*colShifts;
Alinmin = Alinmin+linearBoundsShift;
Alinmax = Alinmax+linearBoundsShift;
%------------------------------------------------------------------%
%  Initialize the Lagrange multipliers on the constraints to zero  %
%------------------------------------------------------------------%
Fmul = zeros(numnonlin+numlin+1,1);
Fstate = Fmul;
%------------------------------------------------------------------%
%   Initialize the Lagrange multipliers on the variables to zero   %
%------------------------------------------------------------------%
xmul = zeros(setup.numvars,1);
xstate = xmul;
%------------------------------------------------------------------%
%       Objrow =1 <=======> First row is the objective row         %
%------------------------------------------------------------------%
ObjRow = 1;
%------------------------------------------------------------------%
%  ObjAdd =0 <====> Do not add anything to the objective function  %
%------------------------------------------------------------------%
ObjAdd = 0;
%------------------------------------------------------------------%
%                     Turn on screen output                        %
%------------------------------------------------------------------%
snscreen on
%------------------------------------------------------------------%
%  Set MYSETUP equal to SETUP for global use in the optimization.  %
%------------------------------------------------------------------%
mysetup = setup;
%------------------------------------------------------------------%
%                   Solve the NLP using SNOPT                      %
%------------------------------------------------------------------% 
Flow = [-Inf; clow; Alinmin];
Fupp = [ Inf; cupp; Alinmax];
[x,F,xmul,Fmul,info,xstate,Fstate,ns,ninf,sinf,mincw,miniw,minrw]...
    = snsolve(xguess,xlow,xupp,xmul,xstate,Flow,Fupp,Fmul,Fstate,...
              ObjAdd,ObjRow,AA,iAfun,jAvar,iGfun,jGvar,userfun);
snprint off;
%------------------------------------------------------------------%
%                Unscale the NLP (Decision) Variables              %
%------------------------------------------------------------------%
result.x = (x-setup.column_shifts)./setup.column_scales;
%------------------------------------------------------------------%
%              Unscale the multipliers on the variables            %
%------------------------------------------------------------------%
result.xmul = setup.Dx*xmul;
%------------------------------------------------------------------%
%            Unscale the multipliers on the constraints            %
%------------------------------------------------------------------%
result.Fmul = setup.DF*Fmul(2:numnonlin+1);
%------------------------------------------------------------------%
%         Extract the multipliers on the linear constraints        %
%------------------------------------------------------------------%
result.Amul = Fmul(numnonlin+2:end);
result.Fmuli = Fmul(1);
result.xstate = xstate;
result.Fstatei = Fstate(1);
result.Fstate = Fstate;
result.Astate = Fstate(numnonlin+2:end);
setup.result = result;
%------------------------------------------------------------------%
%               Record the SNOPT info result                       %
%------------------------------------------------------------------%
setup.SNOPT_info = info;
%------------------------------------------------------------------%
%         Untranscribe the NLP (i.e., convert the optimal          %
%         control problem back to optimal control format)          %
%------------------------------------------------------------------%
setup = gpopsNlp2oc(setup);
