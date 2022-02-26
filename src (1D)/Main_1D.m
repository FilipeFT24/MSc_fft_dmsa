%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format short; beep off; warning('off');
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools_1D.Set_Directories();
%  > Run...
run_1 = true;
% >> ----------------------------------------------------------------------
if run_1
    [inp,msh] = A_1D.WrapUp_A_1D(0.01);
    [msh,pde] = B_1D.WrapUp_B_1D(inp,msh);
end
% >> ----------------------------------------------------------------------

