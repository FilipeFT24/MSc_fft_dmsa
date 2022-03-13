classdef Tools_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------   
        function [] = Set_Directories()
            addpath(genpath('A_1D'));
            addpath(genpath('B_1D'));
            addpath(genpath('C_1D'));
            addpath(genpath('../[Tools]/[Tools - Data]'));
            addpath(genpath('../[Tools]/[Tools - Numerical]'));
            addpath(genpath('../[Tools]/[Tools - Post-processing]'));           
        end
        
        %% > 2. -----------------------------------------------------------      
        function [msh] = Set_msh(msh,stl,s)
            s.stl = stl;
            msh.s = s;
            msh   = orderfields(msh  ,{'d','c','f','s'});
            msh.c = orderfields(msh.c,{'NC','Xc','Vc'});
            msh.f = orderfields(msh.f,{'NF','Xv'});
            msh.s = orderfields(msh.s,{'c','f','bnd_i','bnd_v','stl','A','B','Ac','Bc','xf','xt','Inv'});
        end
    end
end


