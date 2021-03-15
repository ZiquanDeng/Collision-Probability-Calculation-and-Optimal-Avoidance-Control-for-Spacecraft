function y = abs(x)

% AD implementation of abs.m
% Code written by David Benson
% August 2011

y.value = abs(x.value);
outerDerivative = ones(size(x.value(:)));
outerDerivative(x.value(:)<0) = -1;
y = compositeDerivative(x,y,outerDerivative);
y = class(y,'ad');
