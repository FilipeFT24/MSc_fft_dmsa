%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format default;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools.Set_Directory('A');
Tools.Set_Directory('B');
tic;
%  > Class A.
[inp,msh] = A.WrapUp_A;
%  > Class B.
[pde] = 0;%B.WrapUp_B(inp,msh);
%  > ----------------------------------------------------------------------
toc;
% >> --------------------------------------------------------------------
% >> --------------------------------------------------------------------
%  > Working directories.
Tools.Set_Directory('C');
tic;
%  > Class C.
C.WrapUp_C(inp,msh,pde,'blk',0);
%  > --------------------------------------------------------------------
toc;
% >> --------------------------------------------------------------------