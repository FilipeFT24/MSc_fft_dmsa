classdef A1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set up input variables #1.
        function [inp] = Set_inp_1(h)
            %  > msh.
            %    ├─ Examples:
            %        ├─ Example 1: Uniform.
            %        └─ Example 2: Non-uniform.
            %                    ├─ Random.
            %                    └─ Smooth non-uniform.
            %                        ├─  A: Stretching parameter.
            %                        └─  c: Clustered location.
            inp.m.Uniform      = 1;                          %  > Set uniform grid(?).
            if ~inp.m.Uniform
                inp.m.A        = 10.0;                        %  > Stretching parameter.
                inp.m.c        = 0.5;                         %  > Clustered location.
            end
            inp.m.XLim         = [0,1];                       %  > Grid limits(X).
            inp.m.h            = h;                           %  > Grid size.
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Set up input variables #2.
        function [inp] = Set_inp_2(v)
            %  > Analytic function.
            syms x;
            c                  = v(1);
            i                  = v(2);
            f                  = exp(-i.*((x-c).^2));         %  > f.
            inp.f              = matlabFunction(f);           %  > f handle.
            %  > Boundary treatment.
            inp.b.t            = ["Dirichlet","Dirichlet"];   %  > West(1)/East(2).
            if any(inp.b.t ~= "Dirichlet") && any(inp.b.t ~= "Neumann") && any(inp.b.t ~= "Robin")
                return;
            end
            inp.b.change(1)    = 1;                           %  > Add n extra points for boundary face fit.
            inp.b.change(2)    = 0;                           %  > Add n extra points for boundary face fit (error estimator).
            %  > Coefficients.
            inp.c(1)           =  1;                          %  > Convection.
            inp.c(2)           = -1;                          %  > Diffusion.
            %  > Polynomial fit.
            inp.p              = [1,1];                       %  > Convection(1)/Diffusion(2).
            if any(rem(inp.p,2) == 0)                         %  > Allow only p=1,3,5,7,9,etc.
                return;
            end
            %  > P-adaptation.
            %  > #1.
            inp.p_adapt.allow  = 1;                           %  > Allow p-adaptation(?).
            inp.p_adapt.nc     = 50;                          %  > Maximum number of cycles.
            inp.p_adapt.ec_m   = 1.0E-10;                     %  > Minimum (global) discretization error.
            inp.p_adapt.lambda = 0.85;                        %  > Treshold for face selection based on maximum face truncation error (%).
            if ~(inp.p_adapt.lambda < 1)
                return;
            end
            %  > #2.
            inp.p_adapt.opt(1) = 0;                           %  > Use higher-order solution(?).
            inp.p_adapt.opt(2) = 0;                           %  > Add lower-order (predicted) cell truncation error as source term(?).
            if all(inp.p_adapt.opt)                           %  > Only allow w/ lower-order solution.
                return;
            end
            %  > Truncated terms.
            inp.t_terms.allow  = 0;                           %  > Compute truncated terms' magnitude(?).
            inp.t_terms.n      = [3,3];                       %  > Number of terms: Convection(1)/Diffusion(2).
        end
    end
end