%% > Clear memory, clean screen, close any figure.
clear, clc, close all; warning off; beep off;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Set_Directories_2D;
run = 0;
% >> ----------------------------------------------------------------------
if ~run
    h   = 7.5E-02;
    inp = A1_2D.Set_inp([1,1],[1,0.5,0.5]); %  > c/f.
    msh = A2_2D.Set_msh(h);                 %  > h.
    obj = B3_2D.Run_p  (inp,msh);           %  > inp/msh.
else
    obj = C1_2D.Run_p  (1,"2");
end
% >> ----------------------------------------------------------------------
function [] = Set_Directories_2D()
    addpath(genpath('A_2D'));
    addpath(genpath('B_2D'));
    addpath(genpath('C_2D'));
    addpath(genpath('[Post-processing]'));
    addpath(genpath('../[Other stuff]'));
    addpath(genpath('../[Tools]'));
end
% >> ----------------------------------------------------------------------