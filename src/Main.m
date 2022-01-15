%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format default; beep off;
%% > Select test.

tic;
[inp,msh] = Class_1.WrapUp_1;
[bnd,blk] = Class_2.Compute_ErrorPDE(inp,msh);
toc;

Class_3.WrapUp_3(inp,msh,blk,'bnd',0);




