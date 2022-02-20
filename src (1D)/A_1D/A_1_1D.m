classdef A_1_1D
    methods (Static)
        function [inp] = Set_inp(h)
            %% > msh.
            % >> 1. Grid parameters.
            %    ├─ Limits: (Xv)_i,f.
            %    ├─ Examples:
            %        ├─ Example 1: Uniform.
            %                    └─ Uniform/non-uniform grid: h.
            %        └─ Example 2: Non-uniform.
            %                    ├─ Example 2.1. Bulk.
            %                                ├─  2.1.1. Domain percentage: 0 < (Nf)_X,Y < 1.
            %                                └─  2.1.2. Domain stretching: 1 < (Ks)_X,Y < Infinity: e.g.: Ks ~= 3,4,...
            %                    └─ Example 2.2. Wall.
            %                                ├─  2.2.1. Domain percentage: 0 < (Nf)_X,Y < 1.
            %                                ├─  2.2.2. Domain stretching: 1 < (Ks)_X,Y < Infinity. -> e.g.: Ks ~= 1.10,1.01,...
            %                                └─  2.3.2. Location         : e/w
            inp.msh.lim.Xv_i  = 0;
            inp.msh.lim.Xv_f  = 1;
            inp.msh.h         = h;
            inp.msh.eg        = "1";
            inp.msh.s_nu.Nf_X = 0.5;
            inp.msh.s_nu.Ks_X = 5.0;
            inp.msh.s_nu.Lc_X = "e";
            
            %% > pr.
            % >> 2. Problem setup.
            %  > 1. Flow conditions (v,g).
            %  > 2. Boundary conditions (Dirichlet/Neumann/Robin).
            inp.pr.v = 0;
            inp.pr.g = 1;
            inp.pr.w = "Dirichlet";
            inp.pr.e = "Dirichlet";
                       
            %% > fr.
            % >> 3. Flux reconstruction method.
            %  > 1. Test p-adaptation.
            %  > 2. Scheme type (convection/diffusion).
            %    ├─ 2.1. Type 1.
            %        ├─  UDS (Upwind  differencing scheme).
            %        ├─  CDS (Central differencing scheme).
            %        └─  DDS (Downind differencing scheme).
            %    └─ 2.2. Type 2.
            %        └─  Number of neighbours to the left/right.
            inp.fr.p_adapt   = true;
            inp.fr.allow_odd = true;
            inp.fr.n         = 1;
            inp.fr.test_ee   = false;
            inp.fr.type_1.v  = "CDS";
            inp.fr.type_1.g  = "CDS";
            inp.fr.type_2.v  = 1;
            inp.fr.type_2.g  = 1;
        end
    end
end