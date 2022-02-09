classdef Tools_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------   
        function [] = Set_Directories()
            addpath(genpath('A_1D'));
            addpath(genpath('B_1D'));
            addpath(genpath('C_1D'));
            addpath(genpath('D_1D'));
            addpath(genpath('../[Tools]/[Tools - Data]'));
            addpath(genpath('../[Tools]/[Tools - Numerical]'));
            addpath(genpath('../[Tools]/[Tools - Post-processing]'));           
        end
        
        %% > 2. -----------------------------------------------------------      
        function [msh] = Sort_struct(msh)
            % >> msh.
            msh   = orderfields(msh  ,{'d','c','f','s'});
            msh.c = orderfields(msh.c,{'NC','Xc','Vol'});
            msh.f = orderfields(msh.f,{'NF','Xv'});
            msh.s = orderfields(msh.s,{'c','f','bnd','xt','Ls','Tf','xf'});
        end
    end
end


