classdef A_1
    methods (Static)
        %% > Wrap-up A_1.
        function [inp] = WrapUp_A_1()
            % >> 1.
            inp = A_1.Set_inp();
        end
        
        %% > 1. -----------------------------------------------------------
        function [inp] = Set_inp()
            %% > msh.
            % >> 1. Grid limits: (Xv,Yv)_i,(Xv,Yv)_f.
            inp.msh.lim.Xv_i = 0;
            inp.msh.lim.Xv_f = 1;
            inp.msh.lim.Yv_i = 0;
            inp.msh.lim.Yv_f = 1;
            
            % >> 2. Grid type.
            %    ├─ Uniform     grid: h.
            %    └─ Non-uniform grid: Nv=[Nv(X),Nv(Y)].
            inp.msh.h     = 0.05;
            inp.msh.Nv(1) = 25;
            inp.msh.Nv(2) = 25;
            
            % >> 3. Grid types:
            %    ├── v
            %          ├─ Example 1: fft_Distmesh2D (https://github.com/ionhandshaker/distmesh/blob/master/distmesh.m).
            %          └─ Example 2: Random.
            %    └── s
            %          ├─ Example 1: Uniform.
            %          └─ Example 2: Non-uniform.
            %                      ├─ Example 2_1: Bulk.
            %                                  ├─  1. Domain percentage: 0 < (Nf)_X,Y < 1.
            %                                  └─  2. Domain stretching: 1 < (Ks)_X,Y < Infinity: e.g.: Ks ~= 3,4,...
            %                      └─ Example 2_2: Wall.
            %                                  ├─  1. Domain percentage: 0 < (Nf)_X,Y < 1.
            %                                  ├─  2. Domain stretching: 1 < (Ks)_X,Y < Infinity. -> e.g.: Ks ~= 1.10,1.01,...
            %                                  └─  3. Location         : East(E)/West(W), North(N)/South(S).            
            inp.msh.pt       = 'v';
            inp.msh.eg       = '1';
            inp.msh.dm       = '1';
            inp.msh.s_nu.Nf_X = 0.5;
            inp.msh.s_nu.Nf_Y = 0.5;
            inp.msh.s_nu.Ks_X = 7.5;
            inp.msh.s_nu.Ks_Y = 7.5;
            
            %% > pr.
            % >> 4. Flow conditions: 1. Convection parameter: V=[vx,vy].
            %                        2. Diffusion  parameter: G=[gx,gy].
            inp.pr.vx = 0.0;
            inp.pr.vy = 0.0;
            inp.pr.gx = 1.0;
            inp.pr.gy = 1.0;

            %% > fr.
            % >> 5. Flux reconstruction method:
            %  > 1. Simulation type   : 1. Explicit.
            %                           2. Implicit.
            %                           3. Deferred-correction approach (DC).
            %  > 2. Weighting function: 1. Unweighted.
            %                           2. Weighted.
            %  > 3. Face polynomial degree.
            %  > 4. Number of Gauss points/per face.
            %  > 5. Neighbouring type : 1. Vertex (at least 1 common vertex) -> false.
            %                           2. Face   (at least 1 common face  ) -> true.
            %  > 6. Extension type.
            inp.fr.st = 'Implicit';           
            inp.fr.wf = 'Weighted';
            inp.fr.np = 5;
            inp.fr.ng = 5;
            inp.fr.nt = false;
            inp.fr.et = false;
        end
    end
end