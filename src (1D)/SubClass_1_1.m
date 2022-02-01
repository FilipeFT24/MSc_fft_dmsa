classdef SubClass_1_1
    methods (Static)
        %% > Wrap up SubClass_1_1.
        function [msh] = WrapUp_1_1(NC)
            % >> Set up...
            msh.NC = NC;
            msh.NV = NC+1;
            if strcmpi(obj.msh_Type,'Uniform')
                msh = SubClass_1_1.Uniform_mshGenerator(msh);
            elseif strcmpi(obj.msh_Type,'Non-uniform')
                if strcmpi(obj.Stretching_Type,"Random")
                    msh = SubClass_1_1.NonUniform_mshGenerator_1(msh);
                elseif strcmpi(obj.Stretching_Type,"Localized")
                    msh = SubClass_1_1.NonUniform_mshGenerator_2(msh);
                elseif strcmpi(obj.Stretching_Type,"Wall")
                    msh = SubClass_1_1.NonUniform_mshGenerator_3(msh,obj.Wall_Select);
                end
            end
            msh = SubClass_1_1.WrapUp_msh(msh);
        end
        
        % >> #1: Uniform msh.
        function [msh] = Uniform_mshGenerator(msh)
            %  > Xv.
            msh.Xv = linspace(obj.Xv_i,obj.Xv_f,msh.NV);
        end
        
        % >> #2: Non-uniform msh.
        %  > (1): Random msh generator.
        function [msh] = NonUniform_mshGenerator_1(msh)
            %  > Interfaces: i=2,...,NV-1.
            msh.Xv = sort((obj.Xv_f-obj.Xv_i).*rand(msh.NV-2,1)+obj.Xv_i,'ascend');
            %  > Interfaces: i=1,NV.
            msh.Xv = cat(1,obj.Xv_i,msh.Xv,obj.Xv_f); 
        end
        %  > (2): Localized refinement.
        function [msh] = NonUniform_mshGenerator_2(msh)
            % >> Transformation parameters.
            %  > Nf_Unf : Normalized computational domain coordinate (uniform distribution).
            %  > Loc_ptg: Stretching location in "domain percentage".
            %  > B      : Stretching parameter.
            Nf_Unf  = linspace(0,1,msh.NV);
            Loc_ptg = (obj.Nf_Loc-obj.Nf_0)./(1-obj.Nf_0);
            B       = 1./(2.*obj.K_SF).*log((1+(exp(obj.K_SF)-1).*(Loc_ptg))./(1+(exp(-obj.K_SF)-1).*(Loc_ptg)));
            
            %  > Nf: Normalized computational domain coordinate.
            for i = 1:msh.NV
                Nf(i) = (obj.Nf_Loc-obj.Nf_0).*(1+sinh(obj.K_SF.*(Nf_Unf(i)-B))./sinh(obj.K_SF.*B))+obj.Nf_0;
            end
            %  > Xv.
            for i = 1:msh.NV
                msh.Xv(i) = Nf(i).*(obj.Xv_f-obj.Xv_i);
            end
        end
        %  > (3):
        function [msh] = NonUniform_mshGenerator_3(msh,bnd_Loc)
            % >> Transformation parameters.
            %  > Nf_Unf : Normalized computational domain coordinate (uniform distribution).
            Nf_Unf = linspace(0,1,msh.NV);
            Bp     = obj.K_SF+1;
            Bm     = obj.K_SF-1;

            for i = 1:msh.NV
                X  (i) = (Bp./Bm).^(1-Nf_Unf(i));
                Num(i) = Bp-Bm.*X(i);
                Den(i) = X(i)+1;
                x_h(i) = Num(i)./Den(i);
            end
            msh.Xv = obj.Xv_i+(obj.Xv_f-obj.Xv_i).*x_h;
            if strcmpi(bnd_Loc,'West')
                %  > Proceed...
            elseif strcmpi(bnd_Loc,'East')
                %  > Invert mesh.
                for i = 1:msh.NC
                    Vol_inv(i) = msh.Xv(i+1)-msh.Xv(i);
                end
                Vol_inv = flip(Vol_inv);
                for i = 1:msh.NV
                    if i == 1
                        msh.Xv(i) = obj.Xv_i;
                    else
                        msh.Xv(i) = msh.Xv(i-1)+Vol_inv(i-1);
                    end
                end
            end
        end
       
        %% > Tools.
        % >> #0:
        function [] = Add_FolderPaths()
            % >> Add paths...
            addpath [Tools - Data];
            addpath [Tools - Numerical];
            addpath [Tools - Post-processing];
            addpath [Tools - Post-processing]/[Other stuff];
        end
        
        % >> #1:
        function [msh] = WrapUp_msh(msh)
            %  > Xc.
            for i = 1:msh.NC
                msh.Xc(i) = 1./2.*(msh.Xv(i)+msh.Xv(i+1));
            end
            %  > Vol.
            for i = 1:msh.NC
                msh.Vol(i) = msh.Xv(i+1)-msh.Xv(i);
            end
            %  > H_ref.
            msh.H_ref = 1./msh.NC.*(msh.Xv(msh.NV)-msh.Xv(1));
        end       
    end
end