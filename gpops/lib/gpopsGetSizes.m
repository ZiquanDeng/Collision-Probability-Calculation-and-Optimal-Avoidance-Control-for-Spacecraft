function setup = gpopsGetSizes(setup)
%------------------------------------------------------------------%
% Get all sizes in a multiple-phase optimal control problem        %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

%-----------------------------------------------%
% Check field inputs of setup structure         %
%-----------------------------------------------%
if ~isfield(setup,'limits')
    error('"limits" field not found in setup structure')
elseif ~isfield(setup,'name')
    error('"name" field not found in setup structure')
elseif ~isfield(setup,'funcs')
    error('"funcs" field not found in setup structure')
elseif ~isfield(setup,'guess')
    error('"guess" field not found in setup structure')  
end
%-----------------------------------------------%
% Determine the number of phases in the problem %
%-----------------------------------------------%
numphases = length(setup.limits);
limits = setup.limits;
%--------------------------------------------------------------------%
% Get sizes of the following variables in each phase of the problem: %
%--------------------------------------------------------------------%
%    states                                                          %
%    controls                                                        %
%    parameters                                                      %
%    path constraints                                                %
%    event constraints                                               %
%--------------------------------------------------------------------%
sizes = zeros(numphases,5);
for iphase=1:numphases;
    %--------------------------------------------------%
    % Compute the number of nodes in the current phase %
    %--------------------------------------------------%
    setup.limits(iphase).nodes = sum(setup.limits(iphase).nodesPerInterval);
    setup.nodes(iphase) =  setup.limits(iphase).nodes;
    %-----------------------------------------------------%
    % Determine the number of states in the current phase %
    %-----------------------------------------------------%
    if isfield(limits(iphase),'state'),
        if (isfield(limits(iphase).state,'min') && isfield(limits(iphase).state,'max'))
            statesMin = limits(iphase).state.min;
            statesMax = limits(iphase).state.max;
            if (~isempty(statesMin) && ~isempty(statesMax)),
                if isequal(size(statesMin),size(statesMax)),
                    nstates     = size(statesMin,1);
                else
                    errStr = 'State Upper & Lower Bound Matrices Must Be Same Size';
                    error('%s in phase: %i',errStr,iphase)
                end;
                if ~isequal(size(statesMin),[nstates, 3])
                    errStr = 'State Bound Matrices Must Be Number_States x 3';
                    error('%s, in phase: %i',errStr,iphase)
                end
            else
                nstates = 0;
            end;
        else
            errStr = 'Need to Specify Lower and Upper Bounds on States in Phase ';
            error('%s, in phase: %i',errStr,iphase)
        end;
    else
        nstates = 0;
    end;
    %-------------------------------------------------------%
    % Determine the number of controls in the current phase %
    %-------------------------------------------------------%
    if isfield(limits(iphase),'control'),
        if (isfield(limits(iphase).control,'min') && isfield(limits(iphase).control,'max'))
            controlsMin = limits(iphase).control.min;
            controlsMax = limits(iphase).control.max;
            if (~isempty(controlsMin) && ~isempty(controlsMax)),
                if isequal(size(controlsMin),size(controlsMax)),
                    ncontrols   = size(controlsMin,1);
                else
                    errStr = 'Control Upper & Lower Bound Vectors Must Be Same Size';
                    error('%s in phase: %i',errStr,iphase)
                end;
                if ~isequal(size(controlsMin),[ncontrols, 1])
                    errStr = 'Control Bound Vector Must Be Column, Number_Controls x 1';
                    error('%s, in phase: %i',errStr,iphase)
                end
            else
                ncontrols = 0;
            end;
        else
            errStr = 'Need to Specify Lower and Upper Bounds on Controls in Phase ';
            error('%s, in phase: %i',errStr,iphase)
        end;
    else
        ncontrols = 0;
    end;
    %----------------------------------------------------------------%
    % Determine the number of static parameters in the current phase %
    %----------------------------------------------------------------%
    if isfield(limits(iphase),'parameter'),
        if (isfield(limits(iphase).parameter,'min') && isfield(limits(iphase).parameter,'max'))
            parametersMin = limits(iphase).parameter.min;
            parametersMax = limits(iphase).parameter.max;
            if (~isempty(parametersMin) && ~isempty(parametersMax)),
                if isequal(size(parametersMin),size(parametersMax)),
                    nparameters   = size(parametersMin,1);
                else
                    errStr = 'Parameter Upper and Lower Bound Vectors Must Be Same Size';
                    error('%s in phase: %i',errStr,iphase)
                end;
                if ~isequal(size(parametersMin),[nparameters, 1])
                    errStr = 'Parameter Bound Vector Must Be Column, Number_Parameters x 1';
                    error('%s, in phase: %i',errStr,iphase)
                end
            else
                nparameters = 0;
            end;
        else
            errStr = 'Need to Specify Lower and Upper Bounds on Parameters';
            error('%s, in phase: %i',errStr,iphase)
        end;
    else
        nparameters = 0;
    end;
    %-------------------------------------------------%
    % Number of path constraints in the current phase %
    %-------------------------------------------------%
    if isfield(limits(iphase),'path'),
        if (isfield(limits(iphase).path,'min') && isfield(limits(iphase).path,'max'))
            pathsMin = limits(iphase).path.min;
            pathsMax = limits(iphase).path.max;
            if (~isempty(pathsMin) && ~isempty(pathsMax)),
                if isequal(size(pathsMin),size(pathsMax)),
                    npaths        = size(pathsMin,1);
                else
                    errStr = 'Path Upper and Lower Bound Vectors Must Be Same Size';
                    error('%s in phase: %i',errStr,iphase)
                end;
                if ~isequal(size(pathsMin),[npaths, 1])
                    errStr = 'Path Bound Vector Must Be Column, Number_Paths x 1';
                    error('%s, in phase: %i',errStr,iphase)
                end
            else
                npaths = 0;
            end;
        else
            errStr = 'Need to Specify Lower and Upper Bounds on Path Constraints in Phase ';
            error('%s, in phase: %i',errStr,iphase)
        end;
    else
        npaths = 0;
    end;
    %--------------------------------------------------%
    % Number of event constraints in the current phase %
    % -------------------------------------------------%
    if isfield(limits(iphase),'event'),
        if ~isempty(limits(iphase).event),
            if (isfield(limits(iphase).event,'min') && isfield(limits(iphase).event,'max'))
                eventsMin = limits(iphase).event.min;
                eventsMax = limits(iphase).event.max;
                if (~isempty(eventsMin) && ~isempty(eventsMax)),
                    if isequal(size(eventsMin),size(eventsMax)),
                        nevents       = size(eventsMin,1);
                    else
                        errStr = 'Event Upper and Lower Bound Vector Must Be Same Size';
                        error('%s in phase: %i',errStr,iphase)
                    end;
                    if ~isequal(size(eventsMin),[nevents, 1])
                        errStr = 'Event Bound Vector Must Be Column, Number_Events x 1';
                        error('%s, in phase: %i',errStr,iphase)
                    end
                else
                    nevents = 0;
                end;
            else
                errStr = 'Need to Specify Lower and Upper Bounds on Event Constraints in Phase ';
                error('%s, in phase: %i',errStr,iphase)
            end;
        else
            nevents = 0;
        end;
    else
        nevents = 0;
    end;
    %-----------------------------------------------------------------%
    % Set the row vector of sizes in the array SIZES.                 %
    % The IPHASE row of SIZES contains information about phase IPHASE %
    % The order of information in SIZES is as follows:                %
    %   SIZES(IPHASE,1) = number of states in phase IPHASE            %
    %   SIZES(IPHASE,2) = number of controls in phase IPHASE          %
    %   SIZES(IPHASE,3) = number of static parameters in phase IPHASE %
    %   SIZES(IPHASE,4) = number of path constraints in phase  IPHASE %
    %   SIZES(IPHASE,5) = number of event constraints in phase IPHASE %
    % ----------------------------------------------------------------%
    sizes(iphase,:)  = [nstates ncontrols nparameters npaths nevents];
end;
%--------------------------------------------------------------------%
% Get sizes of the linkages in the problem                           %
%--------------------------------------------------------------------%
if ~isfield(setup,'linkages'),
    setup.linkages = [];
end;
linkages = setup.linkages;
numlinkpairs = length(linkages);
numlinks = 0;
for ipair=1:numlinkpairs;
    % ------------------------------------------------------------%
    % Check the sizes of the lower and upper limits               %
    % ------------------------------------------------------------%
    if ~isfield(linkages(ipair),'min') || ~isfield(linkages(ipair),'max'),
        error('Must Specify Both Lower and Upper Bounds on Linkage Constraints');
    else
        linkMin = linkages(ipair).min;
        linkMax = linkages(ipair).max;
        if ~isequal(size(linkMin),size(linkMax))
            error('Linkage Upper and Lower Bound Vector Must Be Same Size');
        end
        if ~isequal(size(linkMin,2), 1),
            errStr = 'Linkage Lower Bound Vector Must Be Column Vector: Number_Links x 1';
            error('%s, in Linkage Pair %i',errStr,ipair)
        end
        if ~isequal(size(linkMax,2), 1),
            errStr = 'Linkage Upper Lower Bound Vector Must Be Column Vector: Number_Links x 1';
            error('%s, in Linkage Pair %i',errStr,ipair)
        end
        numlinks = numlinks + length(linkMin);
    end;
    if ~isfield(linkages(ipair),'left') || ~isfield(linkages(ipair),'right'),
        error('Must Specify Both Phases To be Linked');
    else
        leftPhase = linkages(ipair).left.phase;
        rightPhase = linkages(ipair).right.phase;
        if ~isscalar(leftPhase) || ~isscalar(rightPhase) ...
            || rem(leftPhase,1) ~= 0 || rem(rightPhase,1) ~= 0
            error('Left and Right Phase Numbers Must be Scalar Integers in Linkage Pair %i',ipair)
        end;
    end;
    if isequal(leftPhase,rightPhase) || (leftPhase <= 0) ...
        || (leftPhase > numphases) || (rightPhase <= 0) || (rightPhase > numphases),
        error('Invalid Linkage phase numbers (left: %i, right: %i) in Linkage Pair %i',leftPhase,rightPhase,ipair)
    end
end;
%-------------------------------------------------------%
% Put the array SIZES as a field in th structure SETUP. %
%-------------------------------------------------------%
setup.numphases = numphases;
setup.sizes = sizes;
setup.numlinkpairs = numlinkpairs;
setup.numlinks  = numlinks;

