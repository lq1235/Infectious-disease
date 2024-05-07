clc
clear
a=load('data.txt');
plot(a(:,1),a(:,3)-a(:,4),'s-k','LineWidth',1,'Markersize',13)
hold on
plot(a(:,1),a(:,3)-a(:,4),'k.','LineWidth',1,'Markersize',13)
hold on
% axis([-0.2 3.2,-0.1 15])

xlabel('\fontsize{27} p')
ylabel('\fontsize{27}\DeltaU')
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])