function curvatureRatio = gpopsCurvature(t,x)
%---------------------------------------------------------%
% This function computes the curvature of the trajectory  %
% stored in the variable x with corresponding time points %
% stored in the variable t and outputs the ratio of the   %
% maximum over the mean curvature along the trajectory    %
%---------------------------------------------------------%
tf = t(end);
t0 = t(1);
tau0 = -1;
tauf = +1;
tau = 2*(t-t0)./(tf-t0)+1;
taustep = 2/500;
tstep = (tf-t0)/500;  
tall = [t0:tstep:tf].';
tauall = [-1:taustep:+1].';
tcent = (tall(2:end)+tall(1:end-1))/2;
taucent = (tauall(2:end)+tauall(1:end-1))/2;
% xcent = gpopsLagrange(tcent,t,x);
% xstep = gpopsLagrange(tall,t,x);
xcent = gpopsLagrange(taucent,tau,x);
xstep = gpopsLagrange(tauall,tau,x);
% xderiv = (xcent(2:end,:)-xcent(1:end-1,:))/tstep;
% xdderiv = (xstep(3:end,:)-2*xstep(2:end-1,:)+xstep(1:end-2,:))/tstep^2;
xderiv = (xcent(2:end,:)-xcent(1:end-1,:))/taustep;
xdderiv = (xstep(3:end,:)-2*xstep(2:end-1,:)+xstep(1:end-2,:))/taustep^2;
curvature = abs(xdderiv)./((1+xderiv.^2).^(1.5));
curvatureRatio = max(curvature)/mean(curvature);
