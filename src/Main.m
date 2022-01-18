%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format default;
%% > Run...
%  > Class A.
[inp,msh] = A.WrapUp_A;
%  > Class B.
[pde] = B.WrapUp_B(inp,msh);
%  > Class C.
C.WrapUp_C(inp,msh,pde,'bnd',0);

%  > Measure elapsed time.
% [TA,TB] = Other_stuff.Time_AB(inp,msh);







