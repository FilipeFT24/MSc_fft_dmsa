classdef A_1
    methods (Static)
        %% > Wrap-up A_1.
        % >> ----------------------------------------------------------
        % >> 1. Add folders' path.
        % >> 2. Set 'inp' structure.
        % >> ----------------------------------------------------------
        function [inp] = WrapUp_A_1()
            % >> 1.
            A_1.Add_FolderPaths;
            % >> 2.
            inp = A_1.Set_inp();
        end
        
        %% > 1. -----------------------------------------------------------
        function [] = Add_FolderPaths()
            addpath(genpath('[Tools - Data]'));
            addpath(genpath('[Tools - Numerical]'));
            addpath(genpath('[Tools - Post-processing]'));
        end
        
        %% > 2. -----------------------------------------------------------
        function [inp] = Set_inp()
            %% > msh.
            % >> 1. Grid limits: (Xv,Yv)_i,(Xv,Yv)_f.
            inp.msh.lim.Xv_i = 0;
            inp.msh.lim.Xv_f = 1;
            inp.msh.lim.Yv_i = 0;
            inp.msh.lim.Yv_f = 1;
            
            % >> 2. Vertex coordinates: Nv=[Nv(X),Nv(Y)].
            inp.msh.Nv(1) = 18;
            inp.msh.Nv(2) = 18;
            
            % >> 3. Grid types:
            %  > 1. Type #1.├- v.
            %               └- s.
            %  > 2. Type #2.├- Uniform.
            %               └- Non-uniform.
            %                  ├- Random.
            %                  ├- Bulk.
            %                     ├- 1. Domain percentage: 0 < (Nf)_X,Y < 1.
            %                     └- 2. Domain stretching: 1 < (Ks)_X,Y < Infinity: e.g.: Ks ~= 3,4,...
            %                  └- Wall.
            %                     ├- 1. Domain percentage: 0 < (Nf)_X,Y < 1.
            %                     ├- 2. Domain stretching: 1 < (Ks)_X,Y < Infinity. -> e.g.: Ks ~= 1.10,1.01,...
            %                     └- 3. Location         : East(E)/West(W), North(N)/South(S).
            inp.msh.T_1.t    = 'v';
            inp.msh.T_2.t    = 'Non-uniform';
            inp.msh.T_2.st   = 'Random';
            inp.msh.T_2.Nf_X = 0.5;
            inp.msh.T_2.Nf_Y = 0.5;
            inp.msh.T_2.Ks_X = 5.0;
            inp.msh.T_2.Ks_Y = 5.0;
            
            %% > pr.
            % >> 4. Flow conditions: 1. Convection parameter: V=[vx,vy].
            %                        2. Diffusion  parameter: G=[gx,gy].
            inp.pr.vx = 0.0;
            inp.pr.vy = 0.0;
            inp.pr.gx = 1.0;
            inp.pr.gy = 1.0;
            
            % >> 5. Boundary conditions: 1. EB/EV -> East  boundary type/value.
            %                            2. WB/WV -> West  boundary type/value.
            %                            3. NB/NV -> North boundary type/value.
            %                            4. SB/SV -> South boundary type/value.
            %  > 1. Type.
            inp.pr.t.EB = 'Dirichlet';
            inp.pr.t.WB = 'Dirichlet';
            inp.pr.t.NB = 'Dirichlet';
            inp.pr.t.SB = 'Dirichlet';
            %  > 2. Value (if solution~=analytic).
            inp.pr.v.EB = 0;
            inp.pr.v.WB = 0;
            inp.pr.v.NB = 0;
            inp.pr.v.SB = 0;
            
            %% > fr.
            % >> 6. Flux reconstruction method:
            %  > 1. Simulation type   : 1. Explicit.
            %                           2. Implicit.
            %                           3. Deferred-correction approach (DC).
            %  > 2. Neighbouring type : 1. Vertex (at least 1 common vertex).
            %                           2. Face   (at least 1 common face).
            %  > 3. Stencil extension : 1. Verifiy dimensionless length and extend stencil automatically.
            %                           2. Check stencil limits and extend to the cells within it (Optional).
            %  > 4. Weighting function: 1. Weighted.
            %                           2. Unweighted.
            %  > 5. Face polynomial degree.
            %  > 6. Number of Gauss points/per face.
            %  > 7. Extension type.
            inp.fr.st  = 'Implicit';
            inp.fr.nt  = 'Face';
            inp.fr.ext = 'F';
            inp.fr.wf  = 'Unweighted';
            inp.fr.np  = 9;
            inp.fr.ng  = 1;
            inp.fr.et  = false;
        end
    end
end