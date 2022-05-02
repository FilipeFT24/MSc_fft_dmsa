%% > Clear memory, clean screen, close any figure.
clear, clc, close all; warning off; beep off;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools_1.Set_Directories_2D;
run = [1,0];
% >> ----------------------------------------------------------------------
if run(1)
    inp = A1_2D.Set_inp_2([0.5,0.5,100]); %  > f: (xc,yc)/i.
    msh = A2_2D.Set_msh  (5E-2);          %  > h.
    obj = B3_2D.Run_p    (inp,msh);       %  > inp/msh. 
end
% >> ----------------------------------------------------------------------