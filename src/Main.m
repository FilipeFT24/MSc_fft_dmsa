%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format default; beep off;
%% > Select test.
%  > Add paths...
SubClass_1_1.Add_FolderPaths;

tic;
[inp,msh] = Class_1.Set_Inputs;
[msh,bnd,blk] = Class_2.Compute_ErrorPDE(inp,msh);
toc;

%Fig_0.WrapUp_Fig_0(1,inp.fr.n);
Fig_1.WrapUp_Fig_1(2,inp,msh);
%Fig_2.WrapUp_Fig_2(4,inp,msh,blk);




