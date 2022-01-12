classdef SubClass_1_1
    methods (Static)
        %% > Wrap up SubClass_1_1.
        function [inp,msh] = WrapUp_1_1()
            % >> ----------------------------------------------------------
            % >> 0.   Add folders' path.
            % >> 1.   Set input structure.
            % >> 2.   Uniform mesh.
            % >> 3.   Non-uniform mesh.
            %  > 3.1. Randomly generated mesh (should NOT be used).
            %  > 3.2. Bulk     clustered mesh (based on KS_X, KS_Y, Nf_X, Nf_Y).
            %  > 3.3. Wall     clustered mesh (based on KS_X, KS_Y).
            % >> 4.   Reshape arrays.
            % >> ----------------------------------------------------------
            
            % >> 1.
            inp = SubClass_1_1.Set_inp();
            %  > Select mesh type.
            if strcmpi(inp.msh.T_2.t,'Uniform')
                % >> 2.
                msh = SubClass_1_1.Uniform_mshGenerator(inp);
            elseif strcmpi(inp.msh.T_2.t,'Non-uniform')
                % >> 3.
                if strcmpi(inp.msh.T_2.st,'Random')
                    %  > 3.1.
                    msh = SubClass_1_1.NonUniform_mshGenerator_1(inp);
                elseif strcmpi(inp.msh.T_2.st,'Bulk')
                    %  > 3.2.
                    msh = SubClass_1_1.NonUniform_mshGenerator_2(inp);
                elseif strcmpi(inp.msh.T_2.st,'Wall')
                    %  > 3.3.
                    msh = SubClass_1_1.NonUniform_mshGenerator_3(inp);
                end
            end
        end
        
        %% > Tools.
        %% > 0.) ----------------------------------------------------------
        function [] = Add_FolderPaths()
            addpath(genpath('[Tools - Data]'));
            addpath(genpath('[Tools - Numerical]'));
            addpath(genpath('[Tools - Post-processing]'));
        end
        
        %% > 1.) ----------------------------------------------------------
        function [inp] = Set_inp()
            %% > msh.
            % >> Mesh parameters.
            %  > Limits: (Xv,Yv)_i, (Xv,Yv)_f.
            [inp.msh.lim.Xv_i,inp.msh.lim.Xv_f,inp.msh.lim.Yv_i,inp.msh.lim.Yv_f] = ...
                deal(0,1,0,1);            
            %  > Number of vertices used to generate mesh: (NX,NY).
            [inp.msh.Nv(1),inp.msh.Nv(2)] = ...
                deal(5,5);            
            %  > Mesh type:
            %  - Type #1 (T_1): '^'       -> Triangles.
            %                   's'       -> Squares.
            %  - Type #2 (T_2): 'Uniform' -> Uniform     distribution.
            %                             -> Non-uniform distribution.
            %                                |
            %                                |--> Stretching type (st): --> 1.) 'Random'.
            %                                                               2.) 'Bulk'  .|--> 2.1.) Domain percentage: 0 < (Nf)_X,Y < 1.
            %                                                                                 2.2.) Domain stretching: 1 < (Ks)_X,Y < Infinity. -> e.g.: Ks ~= 3.00,4.00,... 
            %                                                               3.) 'Wall'  .|--> 3.1.) Domain percentage: 0 < (Nf)_X,Y < 1.
            %                                                                                 3.2.) Domain stretching: 1 < (Ks)_X,Y < Infinity. -> e.g.: Ks ~= 1.10,1.01,...
            %                                                                                 3.3.) Location         : East(E)/West(W), North(N)/South(S).
            [inp.msh.T_1.t,inp.msh.T_2.t] = ...
                deal('s','Uniform');
            [inp.msh.T_2.Nf_X,inp.msh.T_2.Nf_Y,inp.msh.T_2.Ks_X,inp.msh.T_2.Ks_Y,inp.msh.T_2.st] = ...
                deal(0.5,0.5,5,5,'Random');

            %% > pr.
            % >> Problem setup.
            %  > Flow conditions    : V(vx,vy) -> Convection parameter.
            %                         G(gx,gy) -> Diffusion  parameter.
            [inp.pr.vx,inp.pr.vy,inp.pr.gx,inp.pr.gy] = ...
                deal(0.1,0,1,0);
            %  > Boundary conditions: EB/EV -> East  boundary type/value.
            %                         WB/WV -> West  boundary type/value.
            %                         NB/NV -> North boundary type/value.
            %                         SB/SV -> South boundary type/value.
            [inp.pr.t.EB,inp.pr.t.WB,inp.pr.t.NB,inp.pr.t.SB] = ...
                deal('Dirichlet','Dirichlet','Dirichlet','Dirichlet');
            [inp.pr.v.EB,inp.pr.v.WB,inp.pr.v.NB,inp.pr.v.SB] = ...
                deal(0,0,0,0); 
            
            %% > fr.
            % >> Flux reconstruction method.
            %  > Integration method: st -> Simulation   type : 1.) Explicit.
            %                                                  2.) Implicit.
            %                                                  3.) Deferred-correction approach.
            %                        nt -> Neighbouring type : 1.) Vertex (at least 1 common vertex).
            %                                                  2.) Face   (at least 1 common face).
            %                        wf -> Weighting function: 1.) Weighted.
            %                                                  2.) Unweighted.
            %                        n  -> Method's order.
            %                        ng -> Number of Gauss points/per face.           
            [inp.fr.st,inp.fr.nt,inp.fr.wf,inp.fr.n,inp.fr.ng] = ...
                deal('Implicit','Face','Unweighted',2,5); 
        end
        
        %% > 2.) ----------------------------------------------------------
        function [msh] = Uniform_mshGenerator(inp)
            % >> Local variables.
            NX_v = inp.msh.Nv(1);
            NY_v = inp.msh.Nv(2);
            Xv_i = inp.msh.lim.Xv_i;
            Xv_f = inp.msh.lim.Xv_f;
            Yv_i = inp.msh.lim.Yv_i;
            Yv_f = inp.msh.lim.Yv_f;
            
            % >> Domain vertices.
            %  > (Xd,Yd).
            Xd_x = linspace(Xv_i,Xv_f,NX_v);
            Yd_y = linspace(Yv_i,Yv_f,NY_v);
            %  > Generate grid.
            [Xd,Yd] = meshgrid(Xd_x,Yd_y);
            %  > (Xv,Yv).
            msh.d.xy_v = SubClass_1_1.Reshape_Arrays(Xd,Yd);
        end
        
        %% > 3.) ----------------------------------------------------------
        % >> 3.1.) --------------------------------------------------------
        function [msh] = NonUniform_mshGenerator_1(inp)
            % >> Local variables.
            NX_v = inp.msh.Nv(1);
            NY_v = inp.msh.Nv(2);
            Xv_i = inp.msh.lim.Xv_i;
            Xv_f = inp.msh.lim.Xv_f;
            Yv_i = inp.msh.lim.Yv_i;
            Yv_f = inp.msh.lim.Yv_f;
            T1   = inp.msh.T_1.t;
                       
            % >> Domain vertices.
            if strcmpi(T1,'s')
                %  > (Xd,Yd).
                Xd_p = sort((Xv_f-Xv_i).*rand(1,NX_v-2)+Xv_i,'ascend');
                Xd_p = cat(2,Xv_i,Xd_p,Xv_f);
                Yd_p = sort((Yv_f-Yv_i).*rand(1,NY_v-2)+Yv_i,'ascend');
                Yd_p = cat(2,Yv_i,Yd_p,Yv_f);
                %  > Generate grid.
                [Xd,Yd] = meshgrid(Xd_p,Yd_p);
            elseif strcmpi(T1,'^')
                %  > (Xd,Yd).
                Xd = sort((Xv_f-Xv_i).*rand(NY_v,NX_v-2)+Xv_i,2,'ascend');
                Xd = cat(2,ones(NY_v,1).*Xv_i,Xd,ones(NY_v,1).*Xv_f);
                Yd = sort((Yv_f-Yv_i).*rand(NY_v-2,NX_v)+Yv_i,1,'ascend');
                Yd = cat(1,ones(1,NX_v).*Yv_i,Yd,ones(1,NX_v).*Yv_f);
            end
            %  > (Xv,Yv).
            msh.d.xy_v = SubClass_1_1.Reshape_Arrays(Xd,Yd);
        end
        % >> 3.2.) --------------------------------------------------------
        function [msh] = NonUniform_mshGenerator_2(inp)
            % >> Local variables.
            %  > Nf_Unf: Normalized computational domain coordinate (uniform distribution).
            %  > Pt    : Stretching location in "domain percentage".
            %  > B     : Stretching parameter.
            NX_v = inp.msh.Nv(1);
            NY_v = inp.msh.Nv(2);
            Xv_i = inp.msh.lim.Xv_i;
            Xv_f = inp.msh.lim.Xv_f;
            Yv_i = inp.msh.lim.Yv_i;
            Yv_f = inp.msh.lim.Yv_f;
            Nf_X = inp.msh.T_2.Nf_X; 
            Nf_Y = inp.msh.T_2.Nf_Y;  
            Ks_X = inp.msh.T_2.Ks_X; 
            Ks_Y = inp.msh.T_2.Ks_Y;  
            
            % >> X.
            Nf_Unf_X = linspace(0,1,NX_v);
            Pt_X     = (Nf_X-0)./(1-0);
            B_X      = 1./(2.*Ks_X).*log((1+(exp(Ks_X)-1).*(Pt_X))./(1+(exp(-Ks_X)-1).*(Pt_X)));
            %  > NF_X.
            i        = 1:NX_v;
            NF_X     = (Nf_X-0).*(1+sinh(Ks_X.*(Nf_Unf_X(i)-B_X))./sinh(Ks_X.*B_X))+0;
            % >> Y.
            Nf_Unf_Y = linspace(0,1,NY_v);
            Pt_Y     = (Nf_Y-0)./(1-0);
            B_Y      = 1./(2.*Ks_Y).*log((1+(exp(Ks_Y)-1).*(Pt_Y))./(1+(exp(-Ks_Y)-1).*(Pt_Y)));
            %  > NF_Y.
            j        = 1:NY_v;
            NF_Y     = (Nf_Y-0).*(1+sinh(Ks_Y.*(Nf_Unf_Y(j)-B_Y))./sinh(Ks_Y.*B_Y))+0;
                       
            % >> Domain vertices.
            %  > (Xd,Yd).
            k    = 1:NX_v;
            l    = 1:NY_v;
            Xd_p = NF_X(k).*(Xv_f-Xv_i);
            Yd_p = NF_Y(l).*(Yv_f-Yv_i);
            %  > Generate grid.
            [Xd,Yd] = meshgrid(Xd_p,Yd_p);
            %  > (Xv,Yv).
            msh.d.xy_v = SubClass_1_1.Reshape_Arrays(Xd,Yd);
        end
        % >> 3.3.) --------------------------------------------------------
        function [msh] = NonUniform_mshGenerator_3(msh)
        end 
        
        %% > 4.) ----------------------------------------------------------
        function [xy_v] = Reshape_Arrays(Xd,Yd)
            xy_v(:,1) = reshape(Xd,[],1);
            xy_v(:,2) = reshape(Yd,[],1);
        end
    end
end