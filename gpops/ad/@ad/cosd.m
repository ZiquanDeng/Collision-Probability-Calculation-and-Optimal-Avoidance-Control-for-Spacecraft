function y = cosd(x);

% AD implementation of cosd.m
% Code written by Ilyssa Sanders and Anil V. Rao
% January 2009

y.value = cos(x.value*pi/180);
outerDerivative = -pi/180*sin(x.value(:)*pi/180);
y = compositeDerivative(x,y,outerDerivative);
y = class(y,'ad');
