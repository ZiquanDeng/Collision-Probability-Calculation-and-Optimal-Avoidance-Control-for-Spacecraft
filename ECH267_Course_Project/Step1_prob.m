clear;
clc;
syms w real;
syms t
% six orbit elements of the spacecraft
e_p=0;   % Eccentricity
M_P=pi/180*-1.5;   % Mean Anomaly
OMEGA_P=pi/180*0;  % longitude of ascending node
omega_P=pi/180*0; % Argument of Perigee
i_P=pi/180*1;  % Orbital inclination ¨the angle between orbit surface and axes surface£©
a_P=6378.004+500;  % semi-length diameter
%six orbit elements of the space object
e_O=0;   %Eccentricity
M_O=pi/180*-1.509;  %Mean Anomaly
OMEGA_O=pi/180*0;   %longitude of ascending node
omega_O=pi/180*0;   %Argument of Perigee
i_O=pi/180*0;   %Orbital inclination ¨the angle between orbit surface and axes surface£©
a_O=6378.004+500;   %semi-length diameter
%dynamics equation DX=A*X+U
A=[0   0   0   1  0   0;
   0   0   0   0  1   0;
   0   0   0   0  0   1;
   0   0   0   0  0  2*w;
   0 -w^2  0   0  0   0;
   0   0 3*w^2 -2*w 0   0];
Phi=statrans(A);   % state tansfer matrix
[mu,a_O,R_P_O,DOTR_P]=gdgshu(e_p,M_P,OMEGA_P,omega_P,i_P,a_P,e_O,M_O,OMEGA_O,omega_O,i_O,a_O);  %The orbit six elements of the two are converted into state quantities and relative state quantities in the geocentric inertial coordinate system
close all
X0=[R_P_O(1) R_P_O(2) R_P_O(3) DOTR_P(1) DOTR_P(2) DOTR_P(3)].';  %The relative state quantity in the geocentric inertial coordinate system is expressed in the centroid orbital coordinate system
w_1=sqrt(mu/a_O^3);  %The angular velocity at the origin of the centroid orbital coordinate system is solved
Phi_t=subs(Phi,{w},{w_1});  %Replace w in the state transition matrix with angular velocity of origin motion in the centroid orbital coordinate system
X=Phi_t*X0;   %According to the state quantity at the initial moment, the state quantity at each ergodic moment is obtained by using the state transition matrix
R=[X(1) X(2) X(3)].';  %The relative distance between a spacecraft and a space object
%The covariance matrix was measured at the initial moment
C_x0=[(0.05/3)^2  0   0   0   0   0;
       0  (0.05/3)^2  0   0   0   0;
       0   0 (0.05/3)^2   0   0   0;
       0    0   0   (0.002/3)^2   0   0;
       0    0   0   0   (0.002/3)^2   0;
       0    0   0   0   0   (0.002/3)^2];
C_x=Phi_t*C_x0*Phi_t.' ; %According to the initial covariance matrix, the covariance matrix of each ergodic moment is obtained by using the state transition matrix
C_r=diag(diag(C_x(1:3,1:3)));   %Take the 3x3 part of the upper left corner of the covariance matrix, that is, the part related to the position quantity
C_r_det=det(C_r);   
syms a b c   
r=[a*sin(b)*cos(c) a*sin(b)*sin(c) a*cos(b)].';  %The position in the inertial coordinate system of the center of mass is represented by polar coordinates
f_R=1/sqrt((2*pi)^3*C_r_det)*exp(-1/2*(r-R).'*inv(C_r)*(r-R))*a^2*sin(b);  %Construct the probability density function
tic

T=35; %The forecast time
fr=zeros(1,T+1);  
% %Taking 1s as the time interval, the collision probabilities at each moment in the approximation process were calculated
for i=22:1:24
    d=i;
    fR=subs(f_R,t,d);
%     f_r=matlabFunction(fR);
    f_r=fR;
    fr(i+1)=jifen3(f_r,0,1,0,pi,0,2*pi);
end
toc

%Draw a figure of collision probability and time
tt=0:1:T; 
plot(tt,fr,'b-')
xlim([0,35])
grid on
set(gca, 'FontSize', 16)
title('Collision probability before taking strategy')
xlabel('Time/s')
ylabel('Collision probability ')
hold on