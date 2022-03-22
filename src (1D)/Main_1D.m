%% > Clear memory, clean screen, close any figure.
clear, clc, close all; warning off; beep off;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools_1D.Set_Directories;
%  > Run...
run_1 = 0;
run_2 = 1;
% >> ----------------------------------------------------------------------
if run_1
    [inp]     = A_1D.Set_A1;
    [msh]     = A_1D.Set_A2(1.0E-2);
    [obj,msh] = B_1D.Run_p (inp,msh);
end
% >> ----------------------------------------------------------------------
if run_2
    C_1D.Check_p(1);
end
% >> ----------------------------------------------------------------------
