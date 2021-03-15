clear;
clc;
syms w real;
syms t
%six orbit elements of the spacecraft
e_p=0;   
M_P=pi/180*-1.5;  
OMEGA_P=pi/180*0;  
omega_P=pi/180*0;  
i_P=pi/180*1;  
a_P=6378.004+500;  
%six orbit elements of the space object
e_O=0;   
M_O=pi/180*-1.509;  
OMEGA_O=pi/180*0;   
omega_O=pi/180*0;   
i_O=pi/180*0;   
a_O=6378.004+500;   
%dynamics equation DX=A*X+U
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
tic

T=35; 
fr=zeros(1,T+1); 

for i=0:1:T
    d=i;
    fR=subs(f_R,t,d);
    f_r=matlabFunction(fR);
    fr(i+1)=integral3(f_r,0,1,0,pi,0,2*pi);  
end


% % Determine the position after maneuver
[aa,bb]=max(fr); 
bb=bb-1;
R_0=matlabFunction(R);
R0=R_0(bb) ;
save('R0.mat','R0')
C_r_0=matlabFunction(C_r); 
C_R=C_r_0(bb); 
save('C_R.mat','C_R')
C_R_det=det(C_R); 


[Rr,feval,exitflag]=fmincon('fmin',R0.',[],[],[],[],[],[],'nonlinear') ;%Solve the position point with the shortest relative distance from the original position when satisfying the critical value of collision probability
save('Rr.mat','Rr')
Phi_bb=subs(Phi_t,t,bb-10); 
Phi_max=double(Phi_bb(1:3,:)); 
save('Phi_max.mat','Phi_max')
syms xx1 xx2 xx3 xx4 xx5 xx6 
xx=[xx1 xx2 xx3 xx4 xx5 xx6];
F=vpa(Phi_max*xx.'-Rr.');  
Dpdx=double(jacobian(F,xx)); 
save('Dpdx.mat','Dpdx')

tao=0; 
X_0=double(subs(X,t,tao)); %The initial value of the state quantity at the start of the maneuver
save('X_0.mat','X_0')

toc