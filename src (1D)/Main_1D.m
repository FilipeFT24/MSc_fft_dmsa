%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format short; beep off; warning('off');
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools_1D.Set_Directories();
%  > Run...
run_1 = 1;
run_2 = 0;
% >> ----------------------------------------------------------------------
if run_1
    [inp,msh] = A_1D.Set_A(0.01);
    [msh,pde] = B_1D.Set_B(inp,msh);
end
% >> ----------------------------------------------------------------------