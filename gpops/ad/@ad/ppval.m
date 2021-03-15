function v=ppval(pp,xx)

% AD implementation of ppval.m
% Code written by David Benson
% Aug 2011
%
% Inputs and Outputs
%      pp: piecewise polynomial data returned by PCHIP, SPLINE, INTERP1, or MKPP
%      xx: values to return (object of AD class)
%      v: polynomail value (object of AD class)
%  (see documentation for MATLAB function ppval)

%eval values
v.value = ppval(pp,xx.value);

% poly derivative
ppDer = pp;
ppDer.order = pp.order-1;
[N M] = size(pp.coefs);
ppDer.coefs = zeros(N,M-1);
for i = 1:M-1
  ppDer.coefs(:,i) = (M-i)*pp.coefs(:,i);
end

% eval derivatives
outerDerivative = ppval(ppDer,xx.value(:));

v = compositeDerivative(xx,v,outerDerivative);
v = class(v,'ad');