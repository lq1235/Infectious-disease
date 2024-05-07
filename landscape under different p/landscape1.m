clc;
clear;

ms = 200;
xmin = [0;0];
xmax = [0.25;1];
x = linspace(xmin(1),xmax(1),ms);
y = linspace(xmin(2),xmax(2),ms);
[X,Y] = meshgrid(x,y);

sample = 'pp23_p=0.00158.txt';
px = load(sample);
p = reshape(px(:,3),ms,ms);
sum(sum(p))
FPx = reshape(px(:,4),ms,ms);
FPy = reshape(px(:,5),ms,ms);
z = trapz(y,trapz(x,p));
Pi = p/z;      
PP = eq(Pi,0)+Pi;     
P_eps=min(min(PP));  
P = P_eps*eq(Pi,0)+Pi;
U = -log(P);
U_max =6.0;
U=U.*(U<U_max)+U_max.*(U>U_max);
surf(X,Y,U)
shading interp;
colormap([jet(256)]);
xlabel('\fontsize{25} P_{R}');
ylabel('\fontsize{25} P_{SS}');
zlabel('\fontsize{25} U')
view([-41.1,7.954]);
zlim([-9.8 20])
set(gca,'xtick',0:0.1:0.25)
set(gca,'ytick',0:0.5:1)
set(gca,'LineWidth',1.2,'Fontsize',25)
set(gca,'TickDir', 'out', 'TickLength', [0.009 0.01])
set(gca,'XTickLabelRotation',0);
set(gca,'YTickLabelRotation',0);
grid on;
set(gca, 'LineWidth', 2.5); 


