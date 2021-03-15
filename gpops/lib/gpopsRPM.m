function RPM = gpopsRPM(numMeshIntervals,meshPoints,nodesPerInterval);

%------------------------------------------------------------------%
% Compute the Legendre-Gauss-Radau (LGR) points, LGR weights, and  %
% LGR differentiation matrix for use with a multiple-interval Radau%
% pseudospectral method.  This code is divided into twoparts:      %
%  Part 1: Compute the LGR points and weights.                     %
%  Part 2: Compute the Radau pseudospectral differentiation matrix %
%------------------------------------------------------------------%
% GPOPS Copyright (c) Anil V. Rao and David Benson                 %
%------------------------------------------------------------------%

%------------------------------------------------------------%
% Compute the LGR points, weights, differentiation matrix,   %
% integration matrix, and unity matrix in each mesh interval %
%------------------------------------------------------------%
for i = 1:numMeshIntervals
  [tau,w,P]     = gpopsGetLGR(nodesPerInterval(i));
  tspan         = meshPoints(i+1)-meshPoints(i);
  tadd          = meshPoints(i+1)+meshPoints(i);
  sSeg{i}       = tspan*(tau+1)/2+meshPoints(i);
  sAll          = [sSeg{i}; meshPoints(i+1)];
  wscaled{i}    = tspan*w/2;
  [D2,Dd2,Do2]  = gpopsCollocD(sAll);
  Dsect{i}      = D2;
  Ddsect{i}     = Dd2;
  Dosect{i}     = Do2;
  Asect{i}      = inv(D2(:,2:end));
  Bsect{i}      = zeros(size(Dsect{i}));
  Bsect{i}(:,1) = 1;
end

if 0,
for i = 1:numMeshIntervals
  [tau,w,P]     = gpopsGetLGR(nodesPerInterval(i));
  tauAll        = [tau; +1];
  tspan         = meshPoints(i+1)-meshPoints(i);
  tadd          = meshPoints(i+1)+meshPoints(i);
  sSeg{i}       = tspan*(tau+1)/2+meshPoints(i);
  sAll          = [sSeg{i}; meshPoints(i+1)];
  wscaled{i}    = tspan*w/2;
  [D2,Dd2,Do2]  = gpopsCollocD(tauAll);
  Dsect{i}      = D2*2/tspan;
  Ddsect{i}     = Dd2*2/tspan;
  Dosect{i}     = Do2*2/tspan;
  Asect{i}      = inv(D2(:,2:end))*tspan/2;
  Bsect{i}      = zeros(size(Dsect{i}));
  Bsect{i}(:,1) = 1;
end
end;

%-----------------------------------------------------------------%
% Compose the points, weights, and matrices computed in each mesh %
% interval into single quantities defined across the entire phase %
%-----------------------------------------------------------------%
s                                = cell2mat(sSeg(:));
w                                = cell2mat(wscaled(:));
A                                = gpopsCompositeA(Asect);
B                                = gpopsCompositeD(Bsect);
RPM.differentiationMatrix        = gpopsCompositeD(Dsect);
RPM.differentiationMatrixDiag    = gpopsCompositeD(Ddsect);
RPM.differentiationMatrixOffDiag = gpopsCompositeD(Dosect);
RPM.integrationMatrix            = gpopsCompositeA(Asect);
RPM.unityMatrix                  = gpopsCompositeD(Bsect);
RPM.Points                       = s;
RPM.Weights                      = w;
