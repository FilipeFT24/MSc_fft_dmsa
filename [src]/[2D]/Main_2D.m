%% > Clear memory, clean screen, close any figure.
clear, clc, close all; warning off; beep off;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Set_Directories_2D;
run = 0;
% >> ----------------------------------------------------------------------
if ~run
    h   = 1.0E-01;                  %  > h.
    s   = [1,1];                  %  > c/f.
    inp = A1_2D.Set_inp(s);       %  > s.
    msh = A2_2D.Set_msh(h,s(2));  %  > h,s(2).
    obj = B3_2D.Run_p  (inp,msh); %  > inp/msh.
else
    obj = C1_2D.Run_p  (0);
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