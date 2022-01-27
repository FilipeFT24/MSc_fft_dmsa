%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format default;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools.Set_Directories();
% >> ----------------------------------------------------------------------
%  > Class A/B.
[inp,msh] = A.WrapUp_A(0.015);
[pde]     = B.WrapUp_B(inp,msh); 
% >> ----------------------------------------------------------------------
%  > Class C/D.
%C.WrapUp_C(inp,msh,pde,'bnd',0);
%D.WrapUp_D(true,true,'1');
% >> ----------------------------------------------------------------------