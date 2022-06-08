%% > Clear memory, clean screen, close any figure.
clear, clc, close all; warning off; beep off;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Set_Directories_1D;
run = [1,0];
% >> ----------------------------------------------------------------------
%  > Standard tests.
if run(1)
    h   = 1.0E-2;
    t   = 1;
    inp = A1_1D.Set_inp([1,2],[50,0.5]); %  > c/f.
    msh = A2_1D.Set_msh(h,t);            %  > h.
    obj = B3_1D.Run_p  (inp,msh);        %  > inp/msh.
end
%  > Other tests.
if run(2)
    B3_1D.Load_p(1);
end
% >> ----------------------------------------------------------------------
function [] = Set_Directories_1D()
    addpath(genpath('A_1D'));
    addpath(genpath('B_1D'));
    addpath(genpath('[Post-processing]'));
    addpath(genpath('../'));
end
% >> ----------------------------------------------------------------------
