clc
clear
a=load('tau-DeltaC.txt')

figure(1)
plot(a(:,1),a(:,2),'s-k','LineWidth',1,'markersize',10)
hold on
plot(a(:,1),a(:,2),'k.','LineWidth',1,'markersize',10)
hold on
xlabel('p','FontSize',27);
ylabel('\tau','FontSize',27);
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
ax = gca();
ax.YRuler.Exponent = 2;
% plot([0.0013 0.0013],[-50 250],'r--','LineWidth',1)
% hold on
% plot([0.0022 0.0022],[-50 250],'b--','LineWidth',1)
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'xtick',0:0.0005:0.003)

figure(2)
plot(a(:,1),a(:,3),'s-k','LineWidth',1,'markersize',10)
hold on
plot(a(:,1),a(:,3),'k.','LineWidth',1,'markersize',10)
hold on
xlabel('p','FontSize',27);
ylabel('\DeltaCC','FontSize',27);
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
ax = gca();
ax.YRuler.Exponent = 3;
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'xtick',0:0.0005:0.0035)