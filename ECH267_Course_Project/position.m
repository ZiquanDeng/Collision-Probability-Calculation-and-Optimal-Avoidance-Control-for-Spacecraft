
clear;%�������
clc;%����
syms w real;
syms t
%�������Ĺ��6����
a_S=6878.004; 
e_S=0;   
i_S=pi/180*1; 
W_S=0;  
w_S=0; 
M_S=pi/180*-1.5;

%�ռ�����Ĺ��6����
a_P=6878.004;   
e_P=0;   
i_P=0;  
W_P=0;   
w_P=0;  
M_P=pi/180*-1.509; 


u=3.986005e+5; %������������

% ��⺽������ռ������ƫ�����E
E=-pi:0.1:pi;
syms t k 
% % ��⺽������E
y_S=inline('E-e_S*sin(E)-M_S','E','e_S','M_S'); %���忪���շ���
yy_S=vectorize(y_S);  %ʹ��������������������
Y_S=feval(yy_S,E,e_S,M_S); 
clf %���ͼ��
plot(E,Y_S,'r');
hold on 
plot(E,zeros(size(E)),'k');
xlabel('t');ylabel('y(t)');
hold off 
zoom on 
[EE_S,yy_S]=ginput(2); %ѡ������ֵ
zoom off
EE_S
EE_SK=EE_S(1)+EE_S(2);
EE_SK
[EEE_S,y_S,exitflag]=fzero(y_S,EE_SK,[],e_S,M_S);
EEE_S

% % ��⺽������E
y_P=inline('E-e_S*sin(E)-M_S','E','e_S','M_S'); %���忪���շ���
yy_P=vectorize(y_P);  %ʹ��������������������
Y_P=feval(yy_P,E,e_P,M_P); 
clf %���ͼ��
plot(E,Y_P,'r');
hold on 
plot(E,zeros(size(E)),'k');
xlabel('t');ylabel('y(t)');
hold off 
zoom on 
[EE_P,yy_P]=ginput(2); %ѡ������ֵ
zoom off
EE_P
EE_PK=EE_P(1)+EE_P(2);
EE_PK
[EEE_P,y_P,exitflag]=fzero(y_P,EE_PK,[],e_P,M_P);
EEE_P

% %��⺽����P�����ǵ㷽��λʸ������Q����ͨ������λʸ������λ�������ٶ����ڹ���ϵ�±�ʾ
PS_1=cos(W_S)*cos(w_S)-sin(W_S)*sin(w_S)*cos(i_S);
PS_2=sin(W_S)*cos(w_S)+cos(W_S)*sin(w_S)*cos(i_S);
PS_3=sin(w_S)*sin(i_S);
QS_1=-cos(W_S)*sin(w_S)-sin(W_S)*cos(w_S)*cos(i_S);
QS_2=-sin(W_S)*sin(w_S)+cos(W_S)*cos(w_S)*cos(i_S);
QS_3=cos(w_S)*sin(i_S);
P_S=[PS_1;PS_2;PS_3];
Q_S=[QS_1;QS_2;QS_3];
rg_S_s=a_S*(cos(EEE_S)-e_S)*P_S+a_S*sqrt(1-e_S^2)*sin(EEE_S)*Q_S ;
rg_S=norm(rg_S_s);
vg_S_s=(sqrt(u*a_S)/rg_S)*(-sin(EEE_S)*P_S+sqrt(1-e_S^2)*cos(EEE_S)*Q_S) ;

% %��⺽����P�����ǵ㷽��λʸ������Q����ͨ������λʸ������λ�������ٶ����ڹ���ϵ�±�ʾ
PP_1=cos(W_P)*cos(w_P)-sin(W_P)*sin(w_P)*cos(i_P);
PP_2=sin(W_P)*cos(w_P)+cos(W_P)*sin(w_P)*cos(i_P);
PP_3=sin(w_P)*sin(i_P);
QP_1=-cos(W_P)*sin(w_P)-sin(W_P)*cos(w_P)*cos(i_P);
QP_2=-sin(W_P)*sin(w_P)+cos(W_P)*cos(w_P)*cos(i_P);
QP_3=cos(w_P)*sin(i_P);
P_P=[PP_1;PP_2;PP_3];
Q_P=[QP_1;QP_2;QP_3];
rg_P_s=a_P*(cos(EEE_P)-e_P)*P_P+a_P*sqrt(1-e_P^2)*sin(EEE_P)*Q_P;
rg_P=norm(rg_P_s);
vg_P_s=(sqrt(u*a_P)/rg_P)*(-sin(EEE_P)*P_P+sqrt(1-e_P^2)*cos(EEE_P)*Q_P);


%��������ϵ���������ϵת��
C=-rg_P_s/norm(rg_P_s);
B=-cross(rg_P_s,vg_P_s)/norm(cross(rg_P_s,vg_P_s));
A=cross(B,C);
T=[A B C]';
rx_S=T*rg_S_s;%Sλ�������������ϵ�б�ʾ
rx_P=T*rg_P_s;%Pλ�������������ϵ�б�ʾ
vx_S=T*vg_S_s;%S�ٶ������������ϵ�б�ʾ
vx_P=T*vg_P_s;%P�ٶ������������ϵ�б�ʾ

%�����Ծ��������ٶ�
rx_S_P=(rx_S-rx_P) ;        
vx_S=(vx_S-vx_P);
close all

X0=[rx_S_P(1) rx_S_P(2) rx_S_P(3) vx_S(1) vx_S(2) vx_S(3)].';  %�����Ĺ�������ϵ�е����״̬����ʾ�����Ĺ������ϵ��


%���״̬ת�ƾ���K
A=[0   0   0   1  0   0;
   0   0   0   0  1   0;
   0   0   0   0  0   1;
   0   0   0   0  0  2*w;
   0 -w^2  0   0  0   0;
   0   0 3*w^2 -2*w 0   0];
syms s;
[M,N]=size(A); %��ȡ����A������������
I=eye(M); %����MxM΢�ĵ�λ��
B=(s*I-A);
B_n=inv(B); 
K=ilaplace(B_n); %��������˹�任��״̬ת�ƾ���

%���K_t������ø�ʱ��״̬��
wz=sqrt(u/a_P^3);  %������Ĺ������ϵԭ���˶����ٶ�
K_t=subs(K,{w},{wz});  %��״̬ת�ƾ����е�w�����Ĺ������ϵԭ���˶����ٶ��滻
X=K_t*X0;   %���ݳ�ʼʱ��״̬��������״̬ת�ƾ�����ø�����ʱ�̵�״̬��
R=[X(1) X(2) X(3)].';  %�������Ϳռ��������Ծ���
%��ʼʱ�̲���Э�������
C_x0=[(0.05/3)^2  0   0   0   0   0;
       0  (0.05/3)^2  0   0   0   0;
       0   0 (0.05/3)^2   0   0   0;
       0    0   0   (0.002/3)^2   0   0;
       0    0   0   0   (0.002/3)^2   0;
       0    0   0   0   0   (0.002/3)^2];
C_x=K_t*C_x0*K_t.'; %���ݳ�ʼʱ��Э�����������״̬ת�ƾ�����ø�����ʱ�̵�Э�������
C_r=diag(diag(C_x(1:3,1:3)));   %ȡЭ����������Ͻ�3x3�Ĳ��֣�����λ������صĲ���
C_r_det=det(C_r);   %��C_r������ʽ


syms a b c   %�������a b c
r=[a*sin(b)*cos(c) a*sin(b)*sin(c) a*cos(b)].';  %�����Ĺ�������ϵ�е�λ���ü������ʾ
f_R=1/sqrt((2*pi)^3*C_r_det)*exp(-1/2*(r-R).'*inv(C_r)*(r-R))*a^2*sin(b);  %��������ܶȺ���
tic

T=35; %Ԥ��ʱ��
fr=zeros(1,36);  %����1x36����
%��1sΪʱ����������ƽ������и�ʱ�̵���ײ����
for i=0:1:35
    d=i;
    fR=subs(f_R,t,d);
    f_r=matlabFunction(fR);
    fr(i+1)=integral3(f_r,0,1,0,pi,0,2*pi);  %������ײ����
end
toc

%����ײ������ʱ���ϵͼ
tt=0:1:35;
plot(tt,fr,'b-')
xlim([0,35])
grid on
set(gca, 'FontSize', 15)
title('����ǰ�ƽ�������ײ����ʾ��ͼ')
xlabel('ʱ��/s')
ylabel('��ײ����')
hold on

[mm,tm]=max(fr); %��ײ�������ֵ����Ӧʱ��
tm=tm-1; %��ײ�������ʱ��
Rsp=matlabFunction(R);
Rspm=Rsp(tm) %��ײ�������ʱ�̺������Ϳռ�������Ծ���
C_r_s=matlabFunction(C_r); 
C_Rm=C_r_s(tm) %��ײ�������ʱ�̵����Э�������
C_Rmd=det(C_Rm); %��ײ�������ʱ�����Э������������ʽ
save('Rspm.mat','Rspm')
save('C_Rm.mat','C_Rm')
[Rs,feval,exitflag]=fmincon('fmin',Rspm.',[],[],[],[],[],[],'nonlinear'); %fmincon����Լ���������Ż����⣬Rs��Ϊ���������λ��
save('Rs.mat','Rs')


% �������Ź�ܲ�������
K_tm=subs(K_t,t,tm-15); %�ӻ�����λ�õ���ײ�������ʱ��λ�õ�״̬ת�ƾ���
K_max=double(K_tm(1:3,:)); %ȡ״̬ת�ƾ����ϰ�3*6���󣨺��ٶȣ��������ģ�4-11��
save('K_max.mat','K_max')
syms xx1 xx2 xx3 xx4 xx5 xx6
xx=[xx1 xx2 xx3 xx4 xx5 xx6];
Y=vpa(K_max*xx.'-Rs.');  %��4-23��
Ydx=double(jacobian(Y,xx)); %YΪ�����ſɱ�����ʽ��ƫ��
save('Ydx.mat','Ydx')
run=0; %������ʼʱ��
X0=double(subs(X,t,run)); %������ʼʱ��״̬����ֵ
save('X0.mat','X0')