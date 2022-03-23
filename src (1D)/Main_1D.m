%% > Clear memory, clean screen, close any figure.
clear, clc, close all; warning off; beep off;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools_1D.Set_Directories;
run = [0,1];
% >> ----------------------------------------------------------------------
if run(1)
    [inp]     = A_1D.Set_A1;
    [msh]     = A_1D.Set_A2(1.0E-3);
    [obj,msh] = B_1D.Run_p (inp,msh);
end
% >> ----------------------------------------------------------------------
if run(2)
    C_1D.Check_P1(1,1);
    C_1D.Check_P2(0,1);
end
% >> ----------------------------------------------------------------------
