%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format default;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools.Set_Directories();
% >> ----------------------------------------------------------------------
%  > Class A/B.
[inp,msh] = A.WrapUp_A;
[pde]     = B.WrapUp_B(inp,msh); 
% >> ----------------------------------------------------------------------
%  > Class C.
C.WrapUp_C(inp,msh,pde,'bnd',0);

pde.E.EN{1}
pde.E.EN{3}


% >> ----------------------------------------------------------------------