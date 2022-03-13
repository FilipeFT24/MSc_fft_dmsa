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
            %    └─ Convection/diffusion coefficients: v,g.
            %% > pr: problem setup.
            %    └─ Differencing schemes order/type.
            %        ├─  p: Polynomial order.
            %        ├─  s: Differencing scheme type.
            %             ├─  "u": UDS (Upwind   biased).
            %             ├─  "c": CDS (Central  biased).
            %             └─  "d": DDS (Downwind biased).
            %% > pa: p-adaptation.
            %    ├─ "p_adapt": Allow p-adaptation.
            %    ├─ "odd"    : Allow odd order differencing schemes.
            %    ├─ "n"      : Number of cycles.
            %    └─ "ee"     : Use error estimator(s) to perform p-adaptation.
            %% > ee: test error estimators.
            %    ├─ "ee"     : Test error estimators (otherwise, perform "p-standard" or "p-adaptation" routines.
            %    └─ Differencing schemes order/type (low/-order/high-order/... solutions).
            %        ├─  p: Polynomial order.
            %        └─  s: Differencing scheme type.
            %             ├─  "u": UDS (Upwind   biased).
            %             ├─  "c": CDS (Central  biased).
            %             └─  "d": DDS (Downwind biased).
            %% > Variables.
            % >> msh.
            inp.msh.Uniform = 0;                     %  > Set uniform grid(?).
            inp.msh.XLim    = [0,1];                 %  > Grid limits.
            inp.msh.h       = h;                     %  > Analytic function.
            inp.msh.A       = 2.5;                   %  > Stretching parameter.
            inp.msh.c       = 0.5;                   %  > Clustered location.
            % >> pv.
            inp.pv.f        = "2";                   %  > Analytic function.
            inp.pv.b(1)     = "Dirichlet";           %  > Left  boundary condition.
            inp.pv.b(2)     = "Dirichlet";           %  > Right boundary condition.
            inp.pv.v(1)     = 1;                     %  > Convection.
            inp.pv.v(2)     = 1;                     %  > Diffusion.
            % >> ps.
            inp.ps.p(1)     = 1;                     %  > Convection.
            inp.ps.p(2)     = 1;                     %  > Diffusion.
            inp.ps.s(1)     = "c";                   %  > Convection.
            inp.ps.s(2)     = "c";                   %  > Diffusion.
            % >> pa.
            inp.pa.adapt    = 0;                     %  > p-adaptation.
            inp.pa.odd      = 0;                     %  > Use UDS/DDS(?).
            inp.pa.n        = 0;                     %  > Number of adaptation cycles.
            inp.pa.ee       = 0;                     %  > Use error estimators to perform adaptation.
            % >> ee.
            inp.ee.test     = 1;                     %  > Test error estimators(?).
            inp.ee.flag     = [0,0,1];               %  > Test flags.            
            if inp.ee.test && nnz(inp.ee.flag) ~= 1
                inp.ee.test = 0;
            else
                if inp.ee.flag(1)
                    inp.ee.p(1)   = 1;               %  > Convection.
                    inp.ee.p(2)   = 1;               %  > Diffusion.
                    inp.ee.s(1)   = "c";             %  > Convection.
                    inp.ee.s(2)   = "c";             %  > Diffusion.
                elseif inp.ee.flag(2)
                elseif inp.ee.flag(3)
                    inp.ee.p(:,1) = [1,2,3];         %  > Convection.
                    inp.ee.p(:,2) = [1,2,3];         %  > Diffusion.
                    inp.ee.s(:,1) = ["c","c","c"];   %  > Convection.
                    inp.ee.s(:,2) = ["c","c","c"];   %  > Diffusion.
                else
                    return;
                end
            end
        end
    end
end