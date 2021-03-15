%Collision Probability after taking avoidance strategy
syms w real;
syms t
load('aa.mat') 

X0=State(end,1:6)';
mu=3.986005e+5; 
a_O=6378.004+500; 
w_1=sqrt(mu/a_O^3); 
A=[0   0   0   1  0   0;
   0   0   0   0  1   0;
   0   0   0   0  0   1;
   0   0   0   0  0  2*w;
   0 -w^2  0   0  0   0;
   0   0 3*w^2 -2*w 0   0]; 
Phi=statrans(A); 
Phi_t=subs(Phi,{w},{w_1});
X=Phi_t*X0; 
r=X(1)^2+X(2)^2+X(3)^2;
f=matlabFunction(sqrt(r));
R=[X(1) X(2) X(3)]';
Phi_t_c=subs(Phi_t,t,t+10); 
C_x0=[(0.05/3)^2  0   0   0   0   0;
       0  (0.05/3)^2  0   0   0   0;
       0   0 (0.05/3)^2   0   0   0;
       0    0   0   (0.002/3)^2   0   0;
       0    0   0   0   (0.002/3)^2   0;
       0    0   0   0   0   (0.002/3)^2]; 
C_x=Phi_t_c*C_x0*Phi_t_c'; 
C_r=diag(diag(C_x(1:3,1:3))); 
C_r_det=det(C_r); 
syms a b c
r=[a*sin(b)*cos(c) a*sin(b)*sin(c) a*cos(b)]';
f_R=1/sqrt((2*pi)^3*C_r_det)*exp(-1/2*(r-R)'*inv(C_r)*(r-R))*a^2*sin(b); 

tic

for i=0:1:25
    d=i;
    fR=subs(f_R,t,d);
    f_r=matlabFunction(fR);
    fr(i+11)=integral3(f_r,0,1,0,pi,0,2*pi); 
end
toc

tt=0:1:35;
plot(tt,fr,'r')
xlim([0,35])
grid on
set(gca, 'FontSize', 16)
title('Collision Probability after Taking Avoidance Strategy')
xlabel('Time/s')
ylabel('Probability')
hold on