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
            inp.A       = 5.0;                            %  > Stretching parameter.
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
            %    └─ Convection/diffusion coefficients: v,g.
            %% > pr: problem setup.
            %    └─ Differencing schemes order/type.
            %        ├─  p: Polynomial order.
            %        └─  s: Polynomial  type (1,2,...).
            %             ├─  -1: UDS biased (1/left).
            %             ├─   0: CDS biased.
            %             └─  +1: DDS biased (1/right).
            %% > pa: p-adaptation.
            %    ├─ p_adapt: Allow p-adaptation(?).
            %    ├─ odd    : Allow UDS/DDS(?).
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
            inp.pv.b       = ["Dirichlet","Dirichlet"];   %  > BC: West(1)/East(2).
            inp.pv.v       = [1,1];                       %  > Coeffs: Convection(1)/Diffusion(2).
            % >> ps.
            inp.ps.p       = [1,1];                       %  > Polynomial order: Convection(1)/Diffusion(2).
            inp.ps.t       = [0,0];                       %  > Polynomial  type: Convection(1)/Diffusion(2).
            inp.ps.add     = 1;
            % >> pa.
            inp.pa.adapt   = 1;                           %  > Allow p-adaptation(?).
            inp.pa.odd     = 0;                           %  > Allow UDS/DDS(?).
            inp.pa.ns      = 3;                           %  > Number of solutions used to "predict" nodal field.
            inp.pa.comp_av = 1;                           %  > Compare w/ analytic values(?).
            % >> tt.
            inp.pt.tt      = 0;                           %  > Compute truncated terms' magnitude(?).
            inp.pt.nt      = [3,3];                       %  > Numb. of terms: Convection(1)/Diffusion(2).
        end       
    end
end