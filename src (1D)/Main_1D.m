%% > Clear memory, clean screen, close any figure.
clear, clc, close all; warning off; beep off;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools_1D.Set_Directories;
run = [1,0];
% >> ----------------------------------------------------------------------
if run(1)
    inp = A1_1D.Set_inp_2([0.5,150]);   %  > f: c/i.
    msh = A2_1D.Set_msh  (1E-2);        %  > h.
    obj = B3_1D.Run_p    (inp,msh);     %  > inp/msh.
end
% >> ----------------------------------------------------------------------
if run(2)
    C_1D.Check_p;
end
% >> ----------------------------------------------------------------------
