%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format short; warning('off');
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools_1D.Set_Directories();
%  > Run...
run_1 = true;
% >> ----------------------------------------------------------------------
if run_1
    %  > Class A/B.
    [inp,msh] = A_1D.WrapUp_A_1D(0.1);
    [msh,pde] = B_1D.WrapUp_B_1D(inp,msh);
    %  > Class C.
    C_1D.WrapUp_C_1D(msh,pde); 
end
% >> ----------------------------------------------------------------------

