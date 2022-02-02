%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format default; warning('off');
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools_1D.Set_Directories();
%  > Run...
run_1 = true;
run_2 = false;
% >> ----------------------------------------------------------------------
if run_1
    %  > Class A/B.
    [inp,msh] = A_1D.WrapUp_A_1D(0.01);
    [pde]     = B_1D.WrapUp_B_1D(inp,msh);
    %  > Class C.
    C_1D.WrapUp_C_1D(msh,pde);
end
% >> ----------------------------------------------------------------------
% >> ----------------------------------------------------------------------
if run_2
    %  > Class D.
    D.WrapUp_D(true,true,'3');
end
% >> ----------------------------------------------------------------------

