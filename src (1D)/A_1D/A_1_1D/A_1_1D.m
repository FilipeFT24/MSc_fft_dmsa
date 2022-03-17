classdef A_1_1D
    methods (Static)
        function [inp] = Set_inp(h)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% > msh: grid parameters/properties.
            %    ├─ Examples:
            %        ├─ Example 1: Uniform.
            %        └─ Example 2: Smooth non-uniform.
            %                    ├─  A: Stretching parameter.
            %                    └─  c: Clustered location.
            %    └─ Limits: XLim = [(Xv)_i,(Xv)_f].
            %% > pv: problem variables.
            %    ├─ Analytic function (f) to be specified in 'B_1_1D':
            %        ├─  1.
            %        ├─  2.
            %        └─  ...
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
            %    └─ n      : Number of cycles.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% > Input variables.
            % >> msh.
            inp.msh.Uniform = 1;                           %  > Set uniform grid(?).
            inp.msh.XLim    = [0,1];                       %  > Grid limits.
            inp.msh.h       = h;                           %  > Grid size.
            inp.msh.A       = 3.5;                         %  > Stretching parameter.
            inp.msh.c       = 0.5;                         %  > Clustered location.
            % >> pv.
            inp.pv.f        = 2;                           %  > f.
            inp.pv.b        = ["Dirichlet","Dirichlet"];   %  > BC: West(1)/East(2).
            inp.pv.v        = [1,1];                       %  > Coeffs: Convection(1)/Diffusion(2).
            % >> ps.
            inp.ps.p        = [1,1];                       %  > Polynomial order: Convection(1)/Diffusion(2).
            inp.ps.t        = [0,0];                       %  > Polynomial  type: Convection(1)/Diffusion(2).
            % >> pa.
            inp.pa.adapt    = 0;                           %  > Allow p-adaptation(?).
            inp.pa.odd      = 0;                           %  > Allow UDS/DDS(?).
            inp.pa.n        = 0;                           %  > Number of cycles.
            %% > Plotting variables.
            % >> pl.
            inp.pl.all      = 1;                           %  > Plot all(?).
            inp.pl.tt       = 0;                           %  > Plot truncated terms w/ analytic solution(?).
            inp.pl.nt       = [5,5];                       %  > Numb. of terms: Convection(1)/Diffusion(2).
        end
    end
end