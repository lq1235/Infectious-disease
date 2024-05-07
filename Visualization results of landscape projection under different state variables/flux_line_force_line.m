clc;
clear;
% Xm=4.0;
% index1=3;
% index2=4;
dif =[0.0000017;0.0000017];
ms = 200;
xmin = [0;0];
xmax = [0.35;0.16]; %注意这里第一个数是x，第二个数是y轴
x = linspace(xmin(1),xmax(1),ms);
y = linspace(xmin(2),xmax(2),ms);
[X,Y] = meshgrid(x,y);

px = load('pp67_p=0.00158.txt');

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
U_max =5;
U=U.*(U<U_max)+U_max.*(U>U_max);

h=figure(1);

view([-21.1133,43.4827]);

%绘制流线、梯度力线
dx = x(2)-x(1);
dy = y(2)-y(1);
[GPx,GPy] = gradient(P,dx,dy);
Jx = FPx.*P - dif(1)*GPx ;
Jy = FPy.*P - dif(2)*GPy ;
E=Jy.^2+Jx.^2;
JJx=Jx./(sqrt(E)+eps);
JJy=Jy./(sqrt(E)+eps);
Fx= dif(1)*GPx./P;
Fy=dif(2)*GPy./P;
F=Fx.^2+Fy.^2;
FFx=Fx./(sqrt(F)+eps);
FFy=Fy./(sqrt(F)+eps);
mg = 5:35:200;
ng = mg;

surf(X,Y,U)
hold on
% 在流和梯度力平面上添加U的等高线
surfc(X, Y, U);
hold on
shading interp;
colormap([jet(256)]);
hold on
%归一化二维流和梯度力
z_coordinate = -39.8;
quiver3(X(mg, ng), Y(mg, ng), z_coordinate * ones(size(X(mg, ng))), JJx(mg, ng), JJy(mg, ng), zeros(size(JJx(mg, ng))), 0.30, 'color', 'r', 'LineWidth', 1);
quiver3(X(mg, ng), Y(mg, ng), z_coordinate * ones(size(X(mg, ng))), FFx(mg, ng), FFy(mg, ng), zeros(size(FFx(mg, ng))), 0.23, 'color', 'k', 'LineWidth', 1);
% quiver3(X(mg, ng), Y(mg, ng), z_coordinate * ones(size(X(mg, ng))), FFx(mg, ng), FFy(mg, ng), zeros(size(FFx(mg, ng))), 'AutoScale','on','AutoScaleFactor',0.2,, 'color', 'k', 'LineWidth', 1);
hold on

xlabel('\fontsize{25} P_{SR}');
ylabel('\fontsize{25} P_{IR}'); %投影不同变量的名称，可以根据可视化绘制的维度修改
zlabel('\fontsize{25} U')
hold on
zlim([-40, max(U(:))+5])


set(gca,'YTickLabelRotation',0)
set(gca,'LineWidth',1.2,'Fontsize',27)
set(gca,'TickDir', 'in', 'TickLength', [0.009 0.01])


