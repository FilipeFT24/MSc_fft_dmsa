%% > Clear memory, clean screen, close any figure.
clear, clc, close all; warning off; beep off;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Set_Directories_1D;
run = [1,0];
% >> ----------------------------------------------------------------------
if run(1)
    inp = A1_1D.Set_inp_2([0.5,1]); %  > f: c/i.
    msh = A2_1D.Set_msh  (2.0E-2);    %  > h.
    obj = B3_1D.Run_p    (inp,msh);   %  > inp/msh.
end
% >> ----------------------------------------------------------------------
function [] = Set_Directories_1D()
    addpath(genpath('A_1D'));
    addpath(genpath('B_1D'));
    addpath(genpath('C_1D'));
    addpath(genpath('[Post-processing]'));
    addpath(genpath('../[Other stuff]'));
    addpath(genpath('../[Tools]'));
end
% >> ----------------------------------------------------------------------
