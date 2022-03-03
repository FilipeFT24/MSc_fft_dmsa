classdef A_1_1D
    methods (Static)
        function [inp] = Set_inp(h)
            %% > msh: mesh parameters.
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
            inp.msh.s_nu.Ks_X = 5;
            inp.msh.s_nu.Lc_X = "e";
            
            %% > pr: problem setup.
            inp.pr.v  = 1;
            inp.pr.g  = 100;
            inp.pr.w  = "Dirichlet";
            inp.pr.e  = "Dirichlet";
            inp.pr.ft = "1";
                       
            %% > fr: flux reconstruction method.
            %  > 1. Test p-adaptation.
            %  > 2. Scheme type (convection/diffusion).
            %    ├─ 2.1. Type 1.
            %        ├─  UDS (Upwind   differencing scheme).
            %        ├─  CDS (Central  differencing scheme).
            %        └─  DDS (Downwind differencing scheme).
            %    └─ 2.2. Type 2.
            %        └─  Number of neighbours to the left/right.
            %  > #1: Standard/p-adaptative routines.
            inp.fr.p_adapt   = true;
            inp.fr.allow_odd = false;
            inp.fr.n         = 1;
            inp.fr.gf        = 1;
            inp.fr.type_1.v  = "UDS";
            inp.fr.type_1.g  = "CDS";
            inp.fr.type_2.v  = 1;
            inp.fr.type_2.g  = 1;
            %  > #2: Test error estimators.
            inp.fr.test_ee   = 0;
        end
    end
end