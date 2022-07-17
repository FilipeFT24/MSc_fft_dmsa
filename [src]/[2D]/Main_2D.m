%% > Clear memory, clean screen, close any figure.
clear, clc, close all; beep off; warning off; opengl hardware; opengl('save','hardware');
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Set_Directories_2D;
% >> ----------------------------------------------------------------------
h   = 7.5E-02;
inp = A1_2D.Set_inp({[1,0,1],1},{[15,3.*pi],[0.50;0.50]}); %  > ch/fh(t,f),wf.
msh = A2_2D.Set_msh(h);                                    %  > h.
obj = B3_2D.Run_p  (inp,msh);                              %  > inp/msh.

tau_c  = obj(end).e.a.n_abs.t.c
tau_f  = obj(end).e.a.n_abs.t.f
e_c    = obj(end).e.a.n_abs.c
nnz_At = obj(end).m.nnz.At

% >> ----------------------------------------------------------------------
function [] = Set_Directories_2D()
    addpath(genpath('A_2D'));
    addpath(genpath('B_2D'));
    addpath(genpath('[Post-processing]'));
    addpath(genpath('../'));
end
% >> ----------------------------------------------------------------------