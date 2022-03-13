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
            %% > Input variables.
            % >> msh.
            inp.msh.Uniform = 0;                     %  > Set uniform grid(?).
            inp.msh.XLim    = [0,2];                 %  > Grid limits.
            inp.msh.h       = h;                     %  > Grid spacing.
            inp.msh.A       = 5.0;                   %  > Stretching parameter.
            inp.msh.c       = 0.5;                   %  > Clustered location.
            % >> pv.
            inp.pv.f        = "1";                   %  > Analytic function.
            inp.pv.b(1)     = "Dirichlet";           %  > BC: West.
            inp.pv.b(2)     = "Dirichlet";           %  > BC: East.
            inp.pv.v(1)     = 1;                     %  > Coeff: Convection.
            inp.pv.v(2)     = 1;                     %  > Coeff: Diffusion.
            % >> ps.
            inp.ps.p(1)     = 1;                     %  > Scheme order: Convection.
            inp.ps.p(2)     = 1;                     %  > Scheme order: Diffusion.
            inp.ps.s(1)     = "c";                   %  > Scheme type:  Convection.
            inp.ps.s(2)     = "c";                   %  > Scheme type:  Diffusion.
            % >> pa.
            inp.pa.adapt    = 0;                     %  > Allow p-adaptation(?).
            inp.pa.odd      = 0;                     %  > Allow UDS/DDS(?).
            inp.pa.n        = 0;                     %  > Number of cycles.
            inp.pa.ee       = 0;                     %  > Use error estimators to perform adaptation(?).
            % >> ee.
            inp.ee.test     = 1;                     %  > Test error estimators(?).
            inp.ee.flag     = [1,0];                 %  > Test flags.            
            if inp.ee.test && nnz(inp.ee.flag) ~= 1
                inp.ee.test = 0;
            else
                if inp.ee.flag(1)
                    inp.ee.p(:,1) = [1,3];           %  > Scheme order: Convection.
                    inp.ee.p(:,2) = [1,3];           %  > Scheme order: Diffusion.
                    inp.ee.s(:,1) = ["c","c"];       %  > Scheme type:  Convection.
                    inp.ee.s(:,2) = ["c","c"];       %  > Scheme type:  Diffusion.
                elseif inp.ee.flag(2)
                    inp.ee.p(:,1) = [1,2,3];         %  > Scheme order: Convection.
                    inp.ee.p(:,2) = [1,2,3];         %  > Scheme order: Diffusion.
                    inp.ee.s(:,1) = ["c","c","c"];   %  > Scheme type:  Convection.
                    inp.ee.s(:,2) = ["c","c","c"];   %  > Scheme type:  Diffusion.
                else
                    return;
                end
            end
            
            %% > Plotting variables.
            % >> pl.
            inp.pl.all   = 1;                        %  > Plot all(?).
            inp.pl.tt    = 1;                        %  > Plot truncated terms w/ analytic solution(?).
            inp.pl.nt(1) = 5;                        %  > Numb. of terms: Convection.
            inp.pl.nt(2) = 5;                        %  > Numb. of terms: Diffusion.
        end
    end
end