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
            %
            %% > pr: problem setup.
            %    ├─ Analytic function (f) to be specified in 'B_1_1D':
            %        ├─  "1".
            %        └─  "2".
            %    ├─ Boundary conditions (west/east boundary):
            %        ├─  "Dirichlet".
            %        ├─  "Neumann".
            %        └─  "Robin".
            %    ├─ Convection/diffusion coefficients: v,g.
            %    ├─ Differencing schemes order/type.
            %        ├─  p: Polynomial order.
            %        ├─  s: Differencing scheme type.
            %             ├─  "u": UDS (Upwind   biased).
            %             ├─  "c": CDS (Central  biased).
            %             └─  "d": DDS (Downwind biased).
            %% > pa: p-adaptation.
            %    ├─ "p_adapt": Allow p-adaptation.
            %    ├─ "odd"    : Allow odd order differencing schemes.
            %    ├─ "n"      : Number of cycles.
            %    ├─ "ee"     : Use error estimator(s) to perform p-adaptation.
            %% > ee: test error estimators.
            %    ├─ "ee"     : Test error estimators (otherwise, perform "p-standard" or "p-adaptation" routines.
            %    ├─ Differencing schemes order/type (low/-order/high-order/... solutions).
            %        ├─  p: Polynomial order.
            %        ├─  s: Differencing scheme type.
            %             ├─  "u": UDS (Upwind   biased).
            %             ├─  "c": CDS (Central  biased).
            %             └─  "d": DDS (Downwind biased).
            %% > Problem variables.
            % >> msh.
            inp.msh.lim.Xv_i = 0;
            inp.msh.lim.Xv_f = 1;
            inp.msh.h        = h;
            inp.msh.t        = "2";
            inp.msh.Nf_X     = 0.50;
            inp.msh.Ks_X     = 5.0;
            inp.msh.Lc_X     = "e";
            % >> ps.
            inp.ps.f         = "2";
            inp.ps.b(1)      = "Dirichlet";
            inp.ps.b(2)      = "Dirichlet"; 
            inp.ps.v(1)      = 1;
            inp.ps.v(2)      = 1;
            inp.ps.p(1)      = 1;
            inp.ps.p(2)      = 1;
            inp.ps.s(1)      = "c";
            inp.ps.s(2)      = "c";
            % >> pa.
            inp.pa.adapt     = 1;
            inp.pa.odd       = 0;
            inp.pa.n         = 0;
            inp.pa.ee        = 0;
            % >> ee.
            inp.ee.test      = 1;
            inp.ee.flag      = [0,0,1];
            inp.ee.p(:,1)    = [1,2];    
            inp.ee.p(:,2)    = [1,2];      
            inp.ee.s(:,1)    = ["c","u"];  
            inp.ee.s(:,2)    = ["c","u"];  
        end
    end
end