classdef A_1_1D
    methods (Static)
        %% > Wrap-up A_1 (1D).
        function [inp] = WrapUp_A_1_1D(h)
            % >> 1.
            inp = A_1_1D.Set_inp(h);
        end
        
        %% > 1. -----------------------------------------------------------
        function [inp] = Set_inp(h)
            %% > msh.
            % >> 1. Grid limits: (Xv)_i,(Xv)_f.
            inp.msh.lim.Xv_i = 0;
            inp.msh.lim.Xv_f = 1;
            
            % >> 2. Grid type.
            %    └─ Uniform/non-uniform grid: h.
            inp.msh.h  = h;
            
            % >> 3. Grid types:
            %    ├─ Example 1: Uniform.
            %          └─ Example 2: Non-uniform.
            %                      ├─ Example 2_1: Bulk.
            %                                  ├─  1. Domain percentage: 0 < (Nf)_X,Y < 1.
            %                                  └─  2. Domain stretching: 1 < (Ks)_X,Y < Infinity: e.g.: Ks ~= 3,4,...
            %                      └─ Example 2_2: Wall.
            %                                  ├─  1. Domain percentage: 0 < (Nf)_X,Y < 1.
            %                                  ├─  2. Domain stretching: 1 < (Ks)_X,Y < Infinity. -> e.g.: Ks ~= 1.10,1.01,...
            %                                  └─  3. Location         : East(E)/West(W), North(N)/South(S).
            inp.msh.eg       = '1';
            inp.msh.s_nu.Nf_X = 0.5;
            inp.msh.s_nu.Ks_X = 5.0;
            
            %% > pr.
            % >> 4. Problem setup.
            %  > 1. Flow conditions    : 1. Convection parameter : v.
            %                            2. Diffusion  parameter : g.
            %  > 2. Boundary conditions: 1. East(E) boundary type: Dirichlet.
            %                            2. West(W) boundary type: Dirichlet/Neumann/Robin.
            inp.pr.v  = 1.0;
            inp.pr.g  = 1.0;
            inp.pr.w  = 'Dirichlet';
            inp.pr.e  = 'Dirichlet';
                       
            %% > fr.
            % >> 5. Flux reconstruction method.
            %  > 1. Simulation type        : 1. Explicit.
            %                                2. Implicit.
            %  > 2. Source term integration: 1. 1D Quadrature (false).
            %                                2. Analytic      (true ).
            %  > 3. Face polynomial degree.
            %  > 4. Number of Gauss points/per face.
            inp.fr.ft = 'Implicit';
            inp.fr.st = false;
            inp.fr.np = 6;
            inp.fr.ng = 3;
        end
    end
end