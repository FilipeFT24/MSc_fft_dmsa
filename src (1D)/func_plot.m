clear; clc;

syms x;
f = exp(-0.1.*((x-0.5).^2));
g = diff(f);
f = matlabFunction(f);
g = matlabFunction(g);

figure;
hold on;
fplot(f,[0,1],'b');
fplot(g,[0,1],'r');
grid on;