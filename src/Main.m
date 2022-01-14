%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format default; beep off;
%% > Select test.

tic;
[inp,msh] = Class_1.WrapUp_1;
[bnd,blk] = Class_2.Compute_ErrorPDE(inp,msh);
toc;

%Fig_0.WrapUp_Fig_0(1,inp.fr.ng);
Fig_1.WrapUp_Fig_1(2,inp,msh,'bnd',0);
%Fig_2.WrapUp_Fig_2(4,inp,msh,blk);




