classdef Tools_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Set_Directories()
            addpath(genpath('A_1D'));
            addpath(genpath('B_1D'));
            addpath(genpath('C_1D'));
            addpath(genpath('D_1D'));
            addpath(genpath('../[Tools]/[Tools - Data]'));
            addpath(genpath('../[Tools]/[Tools - Mesh generation]'));
            addpath(genpath('../[Tools]/[Tools - Numerical]'));
            addpath(genpath('../[Tools]/[Tools - Other stuff]'));
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [msh] = Order_msh(msh)
            msh   = orderfields(msh  ,{'d','c','f','s'});
            msh.c = orderfields(msh.c,{'NC','Xc','Vc'});
            msh.f = orderfields(msh.f,{'NF','Xv'});
            if ~isfield(msh.s,'nt')
                msh.s = orderfields(msh.s,{'c','f','bnd_i','bnd_v','stl','A','B','Ac','Bc','xf','xt','vg','Inv'});
            else
                msh.s = orderfields(msh.s,{'c','f','bnd_i','bnd_v','stl','A','B','Ac','Bc','xf','xt','vg','Inv','nt'});
            end           
        end
        % >> 2.2. ---------------------------------------------------------
        function [pde_e] = Order_pde_e(pde_e)
            pde_e     = orderfields(pde_e    ,{'a','p','eff'});
            pde_e.a   = orderfields(pde_e.a  ,{'c','f','t'});
            pde_e.a.c = orderfields(pde_e.a.c,{'c','c_abs','n','n_abs'});
            pde_e.a.f = orderfields(pde_e.a.f,{'f','f_abs','n','n_abs'});
            pde_e.a.t = orderfields(pde_e.a.t,{'c','c_abs','f','f_abs','n','n_abs'}); 
            pde_e.p   = orderfields(pde_e.p  ,{'c','t'});
            for i = 1:size(pde_e.p.t,2)
                pde_e.p.c{i} = orderfields(pde_e.p.c{i},{'c','c_abs','n','n_abs'});
                pde_e.p.t{i} = orderfields(pde_e.p.t{i},{'c','c_abs','f','f_abs','n','n_abs','s'});
            end
        end
    end
end