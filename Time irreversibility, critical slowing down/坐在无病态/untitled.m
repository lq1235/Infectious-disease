clc
clear

a = [1.3; 0.8472; 1.33; 0.544; 0.515; 0.06];
Xm = 1.0;

% 找出大于1的数的索引
indices = find(a > 1);

% 将大于1的数替换为 2*Xm - 原数值
a(indices) = 2 * Xm - a(indices);
