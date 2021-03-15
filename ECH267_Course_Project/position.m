
clear;%清除数据
clc;%清屏
syms w real;
syms t
%航天器的轨道6根数
a_S=6878.004; 
e_S=0;   
i_S=pi/180*1; 
W_S=0;  
w_S=0; 
M_S=pi/180*-1.5;

%空间物体的轨道6根数
a_P=6878.004;   
e_P=0;   
i_P=0;  
W_P=0;   
w_P=0;  
M_P=pi/180*-1.509; 


u=3.986005e+5; %地心引力常数

% 求解航天器与空间物体的偏近点角E
E=-pi:0.1:pi;
syms t k 
% % 求解航天器的E
y_S=inline('E-e_S*sin(E)-M_S','E','e_S','M_S'); %定义开普勒方程
yy_S=vectorize(y_S);  %使内联函数适于数组运算
Y_S=feval(yy_S,E,e_S,M_S); 
clf %清除图形
plot(E,Y_S,'r');
hold on 
plot(E,zeros(size(E)),'k');
xlabel('t');ylabel('y(t)');
hold off 
zoom on 
[EE_S,yy_S]=ginput(2); %选择坐标值
zoom off
EE_S
EE_SK=EE_S(1)+EE_S(2);
EE_SK
[EEE_S,y_S,exitflag]=fzero(y_S,EE_SK,[],e_S,M_S);
EEE_S

% % 求解航天器的E
y_P=inline('E-e_S*sin(E)-M_S','E','e_S','M_S'); %定义开普勒方程
yy_P=vectorize(y_P);  %使内联函数适于数组运算
Y_P=feval(yy_P,E,e_P,M_P); 
clf %清除图形
plot(E,Y_P,'r');
hold on 
plot(E,zeros(size(E)),'k');
xlabel('t');ylabel('y(t)');
hold off 
zoom on 
[EE_P,yy_P]=ginput(2); %选择坐标值
zoom off
EE_P
EE_PK=EE_P(1)+EE_P(2);
EE_PK
[EEE_P,y_P,exitflag]=fzero(y_P,EE_PK,[],e_P,M_P);
EEE_P

% %求解航天器P（近星点方向单位矢量），Q（半通径方向单位矢量），位置量，速度量在惯性系下表示
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

% %求解航天器P（近星点方向单位矢量），Q（半通径方向单位矢量），位置量，速度量在惯性系下表示
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


%惯性坐标系到相对坐标系转换
C=-rg_P_s/norm(rg_P_s);
B=-cross(rg_P_s,vg_P_s)/norm(cross(rg_P_s,vg_P_s));
A=cross(B,C);
T=[A B C]';
rx_S=T*rg_S_s;%S位置量在相对坐标系中表示
rx_P=T*rg_P_s;%P位置量在相对坐标系中表示
vx_S=T*vg_S_s;%S速度量在相对坐标系中表示
vx_P=T*vg_P_s;%P速度量在相对坐标系中表示

%求解相对距离和相对速度
rx_S_P=(rx_S-rx_P) ;        
vx_S=(vx_S-vx_P);
close all

X0=[rx_S_P(1) rx_S_P(2) rx_S_P(3) vx_S(1) vx_S(2) vx_S(3)].';  %将地心惯性坐标系中的相对状态量表示在质心轨道坐标系中


%求解状态转移矩阵K
A=[0   0   0   1  0   0;
   0   0   0   0  1   0;
   0   0   0   0  0   1;
   0   0   0   0  0  2*w;
   0 -w^2  0   0  0   0;
   0   0 3*w^2 -2*w 0   0];
syms s;
[M,N]=size(A); %获取矩阵A的行数和列数
I=eye(M); %生成MxM微的单位阵
B=(s*I-A);
B_n=inv(B); 
K=ilaplace(B_n); %逆拉普拉斯变换求状态转移矩阵

%求解K_t进而求得各时刻状态量
wz=sqrt(u/a_P^3);  %求解质心轨道坐标系原点运动角速度
K_t=subs(K,{w},{wz});  %将状态转移矩阵中的w用质心轨道坐标系原点运动角速度替换
X=K_t*X0;   %根据初始时刻状态量，利用状态转移矩阵，求得各历经时刻的状态量
R=[X(1) X(2) X(3)].';  %航天器和空间物体的相对距离
%初始时刻测量协方差矩阵
C_x0=[(0.05/3)^2  0   0   0   0   0;
       0  (0.05/3)^2  0   0   0   0;
       0   0 (0.05/3)^2   0   0   0;
       0    0   0   (0.002/3)^2   0   0;
       0    0   0   0   (0.002/3)^2   0;
       0    0   0   0   0   (0.002/3)^2];
C_x=K_t*C_x0*K_t.'; %根据初始时刻协方差矩阵，利用状态转移矩阵，求得各历经时刻的协方差矩阵
C_r=diag(diag(C_x(1:3,1:3)));   %取协方差矩阵左上角3x3的部分，即与位置量相关的部分
C_r_det=det(C_r);   %求C_r的行列式


syms a b c   %定义变量a b c
r=[a*sin(b)*cos(c) a*sin(b)*sin(c) a*cos(b)].';  %将质心惯性坐标系中的位置用极坐标表示
f_R=1/sqrt((2*pi)^3*C_r_det)*exp(-1/2*(r-R).'*inv(C_r)*(r-R))*a^2*sin(b);  %构造概率密度函数
tic

T=35; %预报时长
fr=zeros(1,36);  %建立1x36矩阵
%以1s为时间间隔，计算逼近过程中各时刻的碰撞概率
for i=0:1:35
    d=i;
    fR=subs(f_R,t,d);
    f_r=matlabFunction(fR);
    fr(i+1)=integral3(f_r,0,1,0,pi,0,2*pi);  %计算碰撞概率
end
toc

%画碰撞概率与时间关系图
tt=0:1:35;
plot(tt,fr,'b-')
xlim([0,35])
grid on
set(gca, 'FontSize', 15)
title('机动前逼近过程碰撞概率示意图')
xlabel('时间/s')
ylabel('碰撞概率')
hold on

[mm,tm]=max(fr); %碰撞概率最大值及相应时刻
tm=tm-1; %碰撞概率最大时刻
Rsp=matlabFunction(R);
Rspm=Rsp(tm) %碰撞概率最大时刻航天器和空间物体相对距离
C_r_s=matlabFunction(C_r); 
C_Rm=C_r_s(tm) %碰撞概率最大时刻的误差协方差矩阵
C_Rmd=det(C_Rm); %碰撞概率最大时刻误差协方差矩阵的行列式
save('Rspm.mat','Rspm')
save('C_Rm.mat','C_Rm')
[Rs,feval,exitflag]=fmincon('fmin',Rspm.',[],[],[],[],[],[],'nonlinear'); %fmincon求有约束条件最优化问题，Rs即为所求机动后位置
save('Rs.mat','Rs')


% 能量最优规避策略所用
K_tm=subs(K_t,t,tm-15); %从机动后位置到碰撞概率最大时刻位置的状态转移矩阵
K_max=double(K_tm(1:3,:)); %取状态转移矩阵上半3*6矩阵（含速度），白论文（4-11）
save('K_max.mat','K_max')
syms xx1 xx2 xx3 xx4 xx5 xx6
xx=[xx1 xx2 xx3 xx4 xx5 xx6];
Y=vpa(K_max*xx.'-Rs.');  %（4-23）
Ydx=double(jacobian(Y,xx)); %Y为矩阵，雅可比行列式求偏导
save('Ydx.mat','Ydx')
run=0; %机动开始时刻
X0=double(subs(X,t,run)); %机动开始时刻状态量初值
save('X0.mat','X0')