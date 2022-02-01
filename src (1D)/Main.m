%% > Clear memory, clean screen, close any figure.
clear, clc, close all; format default; warning('off');
%% > Run...
% >> ----------------------------------------------------------------------
%  > Working directories.
Tools_1D.Set_Directories();
%  > Run...
run_1 = true;
run_2 = false;
% >> ----------------------------------------------------------------------
if run_1
    %  > Class A/B.
    [inp,msh] = A_1D.WrapUp_A_1D(0.1);
    [X,Norm,bnd,blk] = Class_2.Compute_ErrorPDE(inp,msh);
    %[pde]     = B_1D.WrapUp_B_1D(inp,msh);
    %  > Class C/D.
    Fig_2_1D.WrapUp_Fig_2_1D(msh,X,Norm);
    %C.WrapUp_C(inp,msh,pde,'bnd',0);
end
% >> ----------------------------------------------------------------------
% >> ----------------------------------------------------------------------
if run_2
    D.WrapUp_D(true,true,'3');
end
% >> ----------------------------------------------------------------------

