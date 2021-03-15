function [mu,a_O,R_P_O,DOTR_P]=gdgshu(P_1,P_2,P_3,P_4,P_5,P_6,E_1,E_2,E_3,E_4,E_5,E_6)
%%The orbital roots are converted to geocentric inertial and relative coordinates
%%P is maneuvering, O is non-maneuvering, and the orbital coordinate system is set on the non-maneuvering, which is fixed with it, that is, on O%%

format long
%%%%%%%%%%%%%%% orbit elements of the spacecraft
E=-pi:0.1:pi;
e_P=P_1; 
M_P=P_2; 
OMEGA_P=P_3; 
omega_P=P_4;
i_P=P_5; 
a_P=P_6; 
mu=3.986005e+5; %Gravity constant

%%orbit elements of the space object
%%%%%%%%%%%%%%not zero(0) is O%%%%%%%%%%%%%
e_O=E_1; 
M_O=E_2; 
OMEGA_O=E_3; 
omega_O=E_4 ; 
i_O=E_5; 
a_O=E_6; 

%%Solve for the auxiliary quantity E of P
y_P=inline('E-e_P*sin(E)-M_P','E','e_P','M_P'); %Define Kepler's equation
y_char_P=vectorize(y_P);
Y_P=feval(y_char_P,E,e_P,M_P); 
clf 
plot(E,Y_P,'r');
hold on 
plot(E,zeros(size(E)),'k')
%[EE_0,y_0,exitflag]=fzero(y,EE(1),[],M0,e0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [EE_P,yy_P]=ginput(1); 
EE_P = 0;
yy_P = 0; 
[EEE_P,y_0,exitflag]=fzero(y_P,EE_P(1),[],e_P,M_P);%%%%waring sequence of e0 M0

%%Solve for the auxiliary quantity E of O
y_O=inline('E-e_O*sin(E)-M_O','E','e_O','M_O'); %Define Kepler's equation
y_char_O=vectorize(y_O);
Y_O=feval(y_char_O,E,e_O,M_O) ;
clf 
plot(E,Y_O,'r');
hold on 
plot(E,zeros(size(E)),'k');
% [EE_0,y_0,exitflag]=fzero(y_O,EE(1),[],M0,e0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EE_O = 0;
yy_O = 0;
[EEE_O,y_0,exitflag]=fzero(y_O,EE_O(1),[],e_O,M_O);%%%%waring sequence of e0 M0


P11_P=cos(OMEGA_P)*cos(omega_P)-sin(OMEGA_P)*sin(omega_P)*cos(i_P);
P21_P=sin(OMEGA_P)*cos(omega_P)+cos(OMEGA_P)*sin(omega_P)*cos(i_P);
P31_P=sin(omega_P)*sin(i_P);

Q11_P=-cos(OMEGA_P)*sin(omega_P)-sin(OMEGA_P)*cos(omega_P)*cos(i_P);
Q21_P=-sin(OMEGA_P)*sin(omega_P)+cos(OMEGA_P)*cos(omega_P)*cos(i_P);
Q31_P=cos(omega_P)*sin(i_P);
P_P=[P11_P;P21_P;P31_P];
Q_P=[Q11_P;Q21_P;Q31_P];

r_0_P=a_P*(cos(EEE_P)-e_P)*P_P+a_P*sqrt(1-e_P^2)*sin(EEE_P)*Q_P ;
r_P=norm(r_0_P);

dot_r_0_P=(sqrt(mu*a_P)/r_P)*(-sin(EEE_P)*P_P+sqrt(1-e_P^2)*cos(EEE_P)*Q_P) ;

P11_O=cos(OMEGA_O)*cos(omega_O)-sin(OMEGA_O)*sin(omega_O)*cos(i_O);
P21_O=sin(OMEGA_O)*cos(omega_O)+cos(OMEGA_O)*sin(omega_O)*cos(i_O);
P31_O=sin(omega_O)*sin(i_O);

Q11_O=-cos(OMEGA_O)*sin(omega_O)-sin(OMEGA_O)*cos(omega_O)*cos(i_O);
Q21_O=-sin(OMEGA_O)*sin(omega_O)+cos(OMEGA_O)*cos(omega_O)*cos(i_O);
Q31_O=cos(omega_O)*sin(i_O);
P_O=[P11_O;P21_O;P31_O];
Q_O=[Q11_O;Q21_O;Q31_O];

r_0_O=a_O*(cos(EEE_O)-e_O)*P_O+a_O*sqrt(1-e_O^2)*sin(EEE_O)*Q_O;
r_O=norm(r_0_O);

dot_r_0_O=(sqrt(mu*a_O)/r_O)*(-sin(EEE_O)*P_O+sqrt(1-e_O^2)*cos(EEE_O)*Q_O);

%%%%%%%%%%%%transfer absolute axes to relative%%%%%%%%%%
%%%%%%Matrix of Oular%%%%%%%%%%%%%%%%%%%%%%%
OULAR_1=-r_0_O/norm(r_0_O);
OULAR_2=-cross(r_0_O,dot_r_0_O)/norm(cross(r_0_O,dot_r_0_O));
OULAR_3=cross(OULAR_2,OULAR_1);
OULAR=[OULAR_3 OULAR_2 OULAR_1]';
R_P=OULAR*r_0_P;                
R_O=OULAR*r_0_O;               
DOT_R_P=OULAR*dot_r_0_P;   
DOT_R_O=OULAR*dot_r_0_O;  

R_P_O=(R_P-R_O) ;        
DOTR_P=(DOT_R_P-DOT_R_O);