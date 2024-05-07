clc
clear
a=load('MFPT.txt');
%MFPT
h=figure(1)
plot(a(:,1),a(:,2),'s-k','LineWidth',1,'markersize',10)
hold on
plot(a(:,1),a(:,2),'k.','LineWidth',1,'markersize',10)
hold on



xlabel("\fontsize{27} p");
ylabel("\fontsize{27} MFPT");

set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.1])


% axis([0.001 0.0025,-0.00005 0.001000001])
xlim([0.0012 0.0023])
ylim([-543080 1343300])
set(gca,'xtick',0:0.0002:0.0040)



h=figure(2)
plot(a(:,1),1./a(:,2),'s-k','LineWidth',1,'markersize',10)
hold on
plot(a(:,1),1./a(:,2),'k.','LineWidth',1,'markersize',10)
hold on

set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.1])
% axis([0 0.00395,-0.00005 0.001000001])
set(gca,'xtick',0:0.0005:0.0040)
xlim([0.001 0.0025])
ylim([-0.0005 0.0013])
xlabel("\fontsize{27} p");
ylabel("\fontsize{27}\nu");
set(gca,'XTickLabelRotation',0);%46是字体的旋转角度
set(gca,'xtick',0:0.00025:0.003)