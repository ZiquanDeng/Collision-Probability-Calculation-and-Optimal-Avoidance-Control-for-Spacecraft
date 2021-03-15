function f=fmin(Rr)
%The cost fuction for calcuating the optimal position to avoid the space
%object
x1=Rr(1);
x2=Rr(2);
x3=Rr(3);
load('R0.mat'); %Read R0
f=(x1-R0(1))^2+(x2-R0(2))^2+(x3-R0(3))^2; %Try to minimize the distance between the optimal position from the orinal position