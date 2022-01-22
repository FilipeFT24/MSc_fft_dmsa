clear; clc;
x = linspace(0,1);
y = linspace(0,1);
[X,Y] = meshgrid(x,y);
Z = (X-0.5).^2+(Y-0.5).^2;
contourf(X,Y,Z,10); axis equal; 

