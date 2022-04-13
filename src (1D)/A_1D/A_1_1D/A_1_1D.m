classdef A_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [inp] = Set_inp_1(h)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% > msh: grid parameters/properties.
            %    ├─ Examples:
            %        ├─ Example 1: Uniform.
            %        └─ Example 2: Smooth non-uniform.
            %                    ├─  A: Stretching parameter.
            %                    └─  c: Clustered location.
            %    └─ Limits: XLim = [(Xv)_i,(Xv)_f].
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% > Input variables.
            % >> msh.
            inp.Uniform = 1;                              %  > Set uniform grid(?).
            inp.XLim    = [0,1];                          %  > Grid limits.
            inp.h       = h;                              %  > Grid size.
            inp.A       = 10.0;                            %  > Stretching parameter.
            inp.c       = 0.5;                            %  > Clustered location.
        end
        % >> 1.2. ---------------------------------------------------------
        function [inp] = Set_inp_2(v)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% > pv: problem variables.
            %    ├─ Analytic function   (f).
            %    ├─ Boundary conditions (W/E boundary):
            %        ├─  "Dirichlet".
            %        ├─  "Neumann".
            %        └─  "Robin".
            %    ├─ Convection/Diffusion coefficients (v/g).
            %    ├─ Differencing schemes order (p).
            %    └─ Add cells to fit boudnary face polynomial (?).
            %% > pa: p-adaptation.
            %    ├─ p_adapt: Allow p-adaptation(?).
            %    └─ n      : Number of solutions used to "predict" nodal field.
            %% > pt: truncated terms.
            %    ├─ tt: Compute truncated terms' magnitude(?).
            %    └─ nt: Number of terms.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% > Input variables.
            % >> pv.
            syms x;
            c              = v(1);
            i              = v(2);
            f              = exp(-i.*((x-c).^2));         %  > f.
            inp.pv.f       = matlabFunction(f);           %  > f handle.
            inp.pv.b       = ["Dirichlet","Robin"];       %  > BC: West(1)/East(2).
            inp.pv.v       = [1,1];                       %  > Coeffs: Convection(1)/Diffusion(2).
            inp.pv.p       = [1,1];                       %  > Polynomial order: Convection(1)/Diffusion(2).
            if any(rem(inp.pv.p,2) == 0)                  %  > Allow only p=1,3,5,etc.
                return;
            end
            inp.pv.add     = 1;                           %  > Add 1 extra point for boundary face fit.
            % >> pa.
            inp.pa.adapt   = 0;                           %  > Allow p-adaptation(?).
            inp.pa.odd     = 0;                           %  > Allow UDS/DDS(?).
            inp.pa.ns      = 1;                           %  > Number of solutions used to "predict" nodal field (>1 to predict "future" field).
            % >> tt.
            inp.pt.tt      = 0;                           %  > Compute truncated terms' magnitude(?).
            inp.pt.nt      = [3,3];                       %  > Numb. of terms: Convection(1)/Diffusion(2).
        end       
    end
end