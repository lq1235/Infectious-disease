clc
clear


a = load('endemic_BarrierHeightChange.txt');
b = load('disease-free_BarrierHeightChange.txt');


X = categorical({'p','\delta','q','r','w'});
X = reordercats(X,{'p','\delta','q','r','w'});


barWidth = 0.35;


offset = barWidth / 2;


X_numeric = double(X);


bar(X_numeric - offset, b, barWidth, 'FaceColor', [0 0.05 0.59]); % 蓝色柱状图
hold on
bar(X_numeric + offset, a, barWidth, 'FaceColor', [1 0.05 0.05]); % 红色柱状图

hold off


set(gca,'ytick',-5:1.0:5.5)
xticks(X_numeric); 
xticklabels(X);    
ylabel("\fontsize{25} \DeltaBarrier")



set(gca,'LineWidth',1.2,'Fontsize',27.4)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])
ylim([-2.5 4.5])

legend('Disease-Free','Endemic',  'FontSize', 17) 