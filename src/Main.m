%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format default;
%% > Select test.
% >> Set directories.
Other_stuff.Add_ABC();

% >> Run...
[inp,msh] = A.WrapUp_A;
[bnd,blk] = B.WrapUp_B(inp,msh);
C.WrapUp_C(inp,msh,blk,'bnd',0);

% >> Measure elapsed time.
% [TA,TB] = Other_stuff.Time_AB(inp,msh);







