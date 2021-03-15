clear;
clc;
syms w real;
syms t
e_p=0;   
M_P=pi/180*-1.5;   
OMEGA_P=pi/180*0;  
omega_P=pi/180*0; 
i_P=pi/180*1;  
a_P=6378.004+500;  

e_O=0;   
M_O=pi/180*-1.509;  
OMEGA_O=pi/180*0;   
omega_O=pi/180*0;   
i_O=pi/180*0;   
a_O=6378.004+500;   

A=[0   0   0   1  0   0;
   0   0   0   0  1   0;
   0   0   0   0  0   1;
   0   0   0   0  0  2*w;
   0 -w^2  0   0  0   0;
   0   0 3*w^2 -2*w 0   0];
Phi=statrans(A);   
[mu,a_O,R_P_O,DOTR_P]=gdgshu(e_p,M_P,OMEGA_P,omega_P,i_P,a_P,e_O,M_O,OMEGA_O,omega_O,i_O,a_O);  
close all
X0=[R_P_O(1) R_P_O(2) R_P_O(3) DOTR_P(1) DOTR_P(2) DOTR_P(3)].';  
w_1=sqrt(mu/a_O^3); 
Phi_t=subs(Phi,{w},{w_1});  
X=Phi_t*X0;   
R=[X(1) X(2) X(3)].'; 

C_x0=[(0.05/3)^2  0   0   0   0   0;
       0  (0.05/3)^2  0   0   0   0;
       0   0 (0.05/3)^2   0   0   0;
       0    0   0   (0.002/3)^2   0   0;
       0    0   0   0   (0.002/3)^2   0;
       0    0   0   0   0   (0.002/3)^2];
C_x=Phi_t*C_x0*Phi_t.' ; 
C_r=diag(diag(C_x(1:3,1:3)));   
C_r_det=det(C_r);   
syms a b c   
r=[a*sin(b)*cos(c) a*sin(b)*sin(c) a*cos(b)].';  
f_R=1/sqrt((2*pi)^3*C_r_det)*exp(-1/2*(r-R).'*inv(C_r)*(r-R))*a^2*sin(b);  

d=23.81771
fR=subs(f_R,t,d);
f_r=matlabFunction(fR);
fr=integral3(f_r,0,1,0,pi,0,2*pi)  









