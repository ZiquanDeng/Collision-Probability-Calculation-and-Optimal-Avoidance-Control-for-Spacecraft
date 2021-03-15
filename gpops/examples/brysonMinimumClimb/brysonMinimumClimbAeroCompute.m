function [CD,CL,eta]=minimumClimbAeroCompute(Mach)

global CONSTANTS

CDdat  = CONSTANTS.CDdat;
CLdat  = CONSTANTS.CLdat;
etadat = CONSTANTS.etadat;

ii     = find(Mach>=0.8);
jj     = find(Mach<0.8);
mpoly  = Mach(ii);
CD  = zeros(length(Mach),1);
CL = zeros(length(Mach),1);
eta = zeros(length(Mach),1);
if ~isempty(ii)
    CD(ii)  = ppval(CDdat,mpoly);
    CL(ii)  = ppval(CLdat,mpoly);
    eta(ii) = ppval(etadat,mpoly);
end

if ~isempty(jj)
    CD(jj)  = 0.013*ones(length(jj),1);
    CL(jj)  = 3.44*ones(length(jj),1);
    eta(jj) = 0.54*ones(length(jj),1);
end