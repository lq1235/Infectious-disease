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
xlabel('t')
ylabel('x')
legend('P_{S}','P_{R}','FontSize',10)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])


%自相关函数
b=load('autocorrelation.txt')
figure(2)
plot(b(:,1),b(:,2),'k-','LineWidth',1,'markersize',10)
hold on
xlabel('t','FontSize',27);
ylabel('A','FontSize',27);
title('Auto correlation');
% axis([-1 70 ,-0.2 1.1])
set(gca,'xtick',0:500:3000)
set(gca,'LineWidth',1.2,'Fontsize',20)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度


%互相关函数
c=load('crosscorrelation.txt')
figure(3)
plot(c(:,1),c(:,2),'k-','LineWidth',1,'markersize',10)
hold on
plot(c(:,1),c(:,3),'r-','LineWidth',1,'markersize',10)
hold on
xlabel('t','FontSize',27);
ylabel('C','FontSize',27);
% ylabel('\tau','FontSize',27);
str={'Time forword C_{XY}','Time reverse C_{YX}'};
legend(str,'FontSize',10)
title('Cross correlation');
% axis([-1 70 ,-0.13 0.13])
set(gca,'xtick',0:500:3000)
set(gca,'LineWidth',1.2,'Fontsize',20)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度


%互相关函数之差
d=load('crosscorrelation-difference.txt')
figure(4)
plot(d(:,1),d(:,2),'k-','LineWidth',1,'markersize',10)
hold on
xlabel('t','FontSize',27);
ylabel('C_{XY}-C_{YX}','FontSize',27);
% str={'Time forword C_{XY}-C_{YX}'};
% legend(str,'FontSize',10)
title('Difference of cross corrlation');
% axis([-1 70 ,-0.2 0.25])
set(gca,'xtick',0:500:3000)
set(gca,'LineWidth',1.2,'Fontsize',20)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])

set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
