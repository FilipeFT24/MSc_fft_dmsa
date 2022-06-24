%% > Clear memory, clean screen, close any figure.
clear, clc, close all; warning off; beep off;
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Set_Directories_2D;
% >> ----------------------------------------------------------------------
h   = 4.0E-02;
inp = A1_2D.Set_inp({[1,2],1},[0.1,0.5,0.5]); %  > c/f.
msh = A2_2D.Set_msh(h);                       %  > h.
obj = B3_2D.Run_p  (inp,msh);                 %  > inp/msh.


tc  = obj(end).e.a.n_abs.t.c
tf  = obj(end).e.a.n_abs.t.f
ec  = obj(end).e.a.n_abs.c
NNZ = obj(end).m{1}.nnz.At


% >> ----------------------------------------------------------------------
function [] = Set_Directories_2D()
    addpath(genpath('A_2D'));
    addpath(genpath('B_2D'));
    addpath(genpath('[Post-processing]'));
    addpath(genpath('../'));
end
% >> ----------------------------------------------------------------------