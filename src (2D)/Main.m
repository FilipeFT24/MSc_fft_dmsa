%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format default; warning('off');
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools.Set_Directories();
%  > Run...
run_1 = true;
run_2 = false;
% >> ----------------------------------------------------------------------
if run_1
    %  > Class A/B.
    [inp,msh] = A.WrapUp_A(0.1);
    [pde]     = B.WrapUp_B(inp,msh);
    %  > Class C/D.
    C.WrapUp_C(inp,msh,pde,'bnd',0);
end
% >> ----------------------------------------------------------------------
% >> ----------------------------------------------------------------------
if run_2
    D.WrapUp_D(true,true,'3');
end
% >> ----------------------------------------------------------------------

