%% > Clear memory, clean screen, close any figure.
clear, clc, close all; warning off; beep off;
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
    [obj,msh] = B_1D.Run_p(inp,msh);
end
% >> ----------------------------------------------------------------------