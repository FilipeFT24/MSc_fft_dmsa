%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format default;
%% > Select test.
%  > "Update" 'obj' class.
obj = Set_NC(obj);
%  > Add paths...
SubClass_1_1.Add_FolderPaths;


%  > Call Class 1.
msh = Class_1.Set_Inputs(obj.NC(1));
%  > Call Class 2.
[X,Norm,bnd,blk] = Class_2.Compute_ErrorPDE(msh);

