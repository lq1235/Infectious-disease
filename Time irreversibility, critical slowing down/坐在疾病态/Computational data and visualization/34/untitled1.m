clc
clear
a=load('series.txt')

%时间序列
t=0:0.1:2999.9
figure(1)
plot(t,a(:,1),'k-','LineWidth',1,'markersize',10)
hold on
plot(t,a(:,2),'r-','LineWidth',1,'markersize',10)
hold on

%自相关函数
b=load('autocorrelation.txt')
figure(2)
plot(b(:,1),b(:,2),'k-','LineWidth',1,'markersize',10)
hold on

%互相关函数
c=load('crosscorrelation.txt')
figure(3)
plot(c(:,1),c(:,2),'k-','LineWidth',1,'markersize',10)
hold on
plot(c(:,1),c(:,3),'k-','LineWidth',1,'markersize',10)
hold on

%互相关函数
d=load('crosscorrelation-difference.txt')
figure(4)
plot(d(:,1),d(:,2),'k-','LineWidth',1,'markersize',10)
hold on



