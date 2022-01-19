%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format default;
%% > Run...
tic;
%  > ----------------------------------------------------------------------
%  > Class A.
[inp,msh] = A.WrapUp_A;
%  > Class B.
[pde] = B.WrapUp_B(inp,msh);
%  > ----------------------------------------------------------------------
toc;
tic;
%  > ----------------------------------------------------------------------
%  > Class C.
C.WrapUp_C(inp,msh,pde,'bnd',0);
%  > ----------------------------------------------------------------------
toc;