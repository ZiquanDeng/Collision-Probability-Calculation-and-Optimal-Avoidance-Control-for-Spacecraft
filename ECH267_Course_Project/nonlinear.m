function [g,h]=nonlinear(Rr)
%%Nonlinear constraint of fmincon
x1=Rr(1);
x2=Rr(2);
x3=Rr(3);
load('C_R.mat'); %Read the error covariance matrix at the moment of maximum collision probability
syms x y z
r=[x*sin(y)*cos(z) x*sin(y)*sin(z) x*cos(y)].';
R_r=[x1 x2 x3].';
C_R_det=det(C_R); 
f_0=1/sqrt((2*pi)^3*C_R_det)*exp(-1/2*(r-R_r).'*inv(C_R)*(r-R_r))*x^2*sin(y); %Construct the probability density function
f_1=matlabFunction(f_0);
f_2= integral3(f_1,0,1,0,pi,0,2*pi); 
g=f_2-1e-7; %Nonlinear inequality constraints
h=[]; %Nonlinear equality constraint
