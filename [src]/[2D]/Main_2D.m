%% > Clear memory, clean screen, close any figure.
clear, clc, close all; warning off; beep off;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Set_Directories_2D;
% >> ----------------------------------------------------------------------
h   = 5.0E-02;
inp = A1_2D.Set_inp({[1,1],1},[75,0.5,0.5]); %  > c/f.
msh = A2_2D.Set_msh(h);                      %  > h.
obj = B3_2D.Run_p  (inp,msh);                %  > inp/msh.
% >> ----------------------------------------------------------------------
function [] = Set_Directories_2D()
    addpath(genpath('A_2D'));
    addpath(genpath('B_2D'));
    addpath(genpath('C_2D'));
    addpath(genpath('[Post-processing]'));
    addpath(genpath('../'));
end
% >> ----------------------------------------------------------------------