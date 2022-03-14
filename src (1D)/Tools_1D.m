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
        % >> 2.1. ---------------------------------------------------------
        function [msh] = Order_msh(msh)
            msh   = orderfields(msh  ,{'d','c','f','s'});
            msh.c = orderfields(msh.c,{'NC','Xc','Vc'});
            msh.f = orderfields(msh.f,{'NF','Xv'});
            if ~(isfield(msh.s,'p') && isfield(msh.s,'nt'))
                msh.s = orderfields(msh.s,{'c','f','bnd_i','bnd_v','stl','A','B','Ac','Bc','xf','xt','vg','Inv'});
            else
                msh.s = orderfields(msh.s,{'c','f','bnd_i','bnd_v','stl','A','B','Ac','Bc','xf','xt','Inv','vg','p','nt'});
            end           
        end
        % >> 2.2. ---------------------------------------------------------
        function [pde_e] = Order_pde_e(pde_e)
            pde_e   = orderfields(pde_e  ,{'c','f','t'});
            pde_e.c = orderfields(pde_e.c,{'c','c_abs','n','n_abs'});
            pde_e.f = orderfields(pde_e.f,{'f','f_abs','n','n_abs'});
            if ~isfield(pde_e.t,'a')
                pde_e.t = orderfields(pde_e.t,{'c','c_abs','f','f_abs','n','n_abs'});
            else
                pde_e.t = orderfields(pde_e.t,{'a','c','c_abs','f','f_abs','n','n_abs'});
            end
        end
    end
end