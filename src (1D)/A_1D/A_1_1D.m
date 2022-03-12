classdef A_1_1D
    methods (Static)
        function [inp] = Set_inp(h)
            %% > msh: grid parameters/properties.
            %    ├─ Examples:
            %        ├─ Example 1: Uniform.
            %        └─ Example 2: Smoothly non-uniform.
            %                    ├─  A: Stretching parameter.
            %                    └─  c: Clustered location.
            %    └─ Limits: XLim = [(Xv)_i,(Xv)_f].
            %% > pv: problem variables.
            %    ├─ Analytic function (f) to be specified in 'B_1_1D':
            %        ├─  "1".
            %        ├─  "2".
            %        └─  ...
            %    ├─ Boundary conditions (west/east boundary):
            %        ├─  "Dirichlet".
            %        ├─  "Neumann".
            %        └─  "Robin".
            %    ├─ Convection/diffusion coefficients: v,g.
            %% > pr: problem setup.
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
            %% > Variables.
            % >> msh.
            inp.msh.Uniform = 0;
            inp.msh.XLim    = [0,1];
            inp.msh.h       = h;   
            inp.msh.A       = 5.0;
            inp.msh.c       = 0.5;
            % >> pv.
            inp.pv.f        = "1";
            inp.pv.b(1)     = "Dirichlet";
            inp.pv.b(2)     = "Dirichlet"; 
            inp.pv.v(1)     = 1;
            inp.pv.v(2)     = 1;
            % >> ps.
            inp.ps.p(1)     = 1;
            inp.ps.p(2)     = 1;
            inp.ps.s(1)     = "c";
            inp.ps.s(2)     = "c";
            % >> pa.
            inp.pa.adapt    = 0;
            inp.pa.odd      = 0;
            inp.pa.n        = 0;
            inp.pa.ee       = 0;
            % >> ee.
            inp.ee.test     = 0;
            inp.ee.flag     = [0,0,1];
            inp.ee.p(:,1)   = [1,2];    
            inp.ee.p(:,2)   = [1,2];      
            inp.ee.s(:,1)   = ["c","c"];  
            inp.ee.s(:,2)   = ["c","c"];  
        end
    end
end