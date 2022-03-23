%% > Clear memory, clean screen, close any figure.
clear, clc, close all; warning off; beep off;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools_1D.Set_Directories;
run = [1,0];
% >> ----------------------------------------------------------------------
if run(1)
    [inp]     = A_1D.Set_A1;
    [msh]     = A_1D.Set_A2(0.25E-2);
    [obj,msh] = B_1D.Run_p (inp,msh);
end
% >> ----------------------------------------------------------------------
if run(2)
    C_1D.Check_P1(1,1);
end
% >> ----------------------------------------------------------------------
