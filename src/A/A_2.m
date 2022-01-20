classdef A_2
    methods (Static)
        %% > Wrap-up A_2.
        % >> --------------------------------------------------------------
        % >> 1.   Uniform grid.
        % >> 2.   Non-uniform grid.
        %  > 2.1. Grid: Randomly generated point cloud.
        %  > 2.2. Grid: Blust clustering.
        %  > 2.3. Grid: Wall  clustering.
        % >> 3.   Reshape arrays.
        % >> --------------------------------------------------------------
        function [msh] = WrapUp_A_2(inp)
            % >> Local variables.
            T_2    = inp.msh.T_2.t;
            T_2_st = inp.msh.T_2.st;
            
            %  > Select grid type.
            if strcmpi(T_2,'Uniform')
                % >> 1.
                msh = A_2.Uniform_mshGenerator(inp);
            elseif strcmpi(T_2,'Non-uniform')
                % >> 2.
                if strcmpi(T_2_st,'Random')
                    %  > 2.1.
                    msh = A_2.NonUniform_mshGenerator_1(inp);
                elseif strcmpi(T_2_st,'Bulk')
                    %  > 2.2.
                    msh = A_2.NonUniform_mshGenerator_2(inp);
                elseif strcmpi(T_2_st,'Wall')
                end
            end
        end
        
        %% > 1. -----------------------------------------------------------
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
            msh.d.xy_v = A_2.Reshape_Arrays(Xd,Yd);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [msh] = NonUniform_mshGenerator_1(inp)
            % >> Local variables.
            NX_v = inp.msh.Nv(1);
            NY_v = inp.msh.Nv(2);
            Xv_i = inp.msh.lim.Xv_i;
            Xv_f = inp.msh.lim.Xv_f;
            Yv_i = inp.msh.lim.Yv_i;
            Yv_f = inp.msh.lim.Yv_f;
            T1   = inp.msh.T_1.t;
            % rng default;
            
            % >> Domain vertices.
            if strcmpi(T1,'s')
                %  > (Xd,Yd).
                Xd_p = sort((Xv_f-Xv_i).*rand(1,NX_v-2)+Xv_i,'ascend');
                Xd_p = cat(2,Xv_i,Xd_p,Xv_f);
                Yd_p = sort((Yv_f-Yv_i).*rand(1,NY_v-2)+Yv_i,'ascend');
                Yd_p = cat(2,Yv_i,Yd_p,Yv_f);
                %  > Generate grid.
                [Xd,Yd] = meshgrid(Xd_p,Yd_p);
            elseif strcmpi(T1,'v')
                %  > (Xd,Yd).
                Xd = sort((Xv_f-Xv_i).*rand(NY_v,NX_v-2)+Xv_i,2,'ascend');
                Xd = cat(2,ones(NY_v,1).*Xv_i,Xd,ones(NY_v,1).*Xv_f);
                Yd = sort((Yv_f-Yv_i).*rand(NY_v-2,NX_v)+Yv_i,1,'ascend');
                Yd = cat(1,ones(1,NX_v).*Yv_i,Yd,ones(1,NX_v).*Yv_f);
            end
            %  > (Xv,Yv).
            msh.d.xy_v = A_2.Reshape_Arrays(Xd,Yd);
        end
        % >> 2.2. ---------------------------------------------------------
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
            msh.d.xy_v = A_2.Reshape_Arrays(Xd,Yd);
        end
        % >> 2.3. ---------------------------------------------------------
        function [msh] = NonUniform_mshGenerator_3(msh)
        end
        
        %% > 3. -----------------------------------------------------------
        function [xy_v] = Reshape_Arrays(Xd,Yd)
            xy_v(:,1) = reshape(Xd,[],1);
            xy_v(:,2) = reshape(Yd,[],1);
        end
    end
end