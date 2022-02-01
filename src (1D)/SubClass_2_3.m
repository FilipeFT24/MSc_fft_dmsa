classdef SubClass_2_3
    methods (Static)
        %% > Wrap up SubClass_2_3.
        function [msh_o,X_o,Norm_o] = WrapUp_2_3(msh_i,X_i,Norm_i)
           [msh_o,X_o,Norm_o] = SubClass_2_3.Select_EligibleCells(msh_i,X_i,Norm_i);
        end
        
        %% > Step #1: Select eligible cells for refinement.
        function [msh_o,X_o,Norm_o] = Select_EligibleCells(msh_i,X_i,Norm_i)
            %  > Error treshold.
            Er_Trsh = obj.K_ref.*cell2mat(Norm_i.E(1));
            %  > Index of the selected cells.
            msh_o.j_ref = find(X_i.Error > Er_Trsh);
            
            for i = 1:msh_i.NC
                %  > Lhs and rhs cell vertices.
                i_Ini = msh_i.Xv(i);
                i_Fin = msh_i.Xv(i+1);
                %  > Store...
                msh_o.Xv_i{i} = [i_Ini,i_Fin];
                if ismember(i,msh_o.j_ref)
                    msh_o.Xv_i{i} = linspace(i_Ini,i_Fin,3);
                end 
            end
            %  > Xv.
            msh_o.Xv = cell2mat(msh_o.Xv_i);
            msh_o.Xv = unique(msh_o.Xv);
            %  > NV,NC.
            msh_o.NV = length(msh_o.Xv);
            msh_o.NC = msh_o.NV-1;
            %  > Xc,Vol,H.
            msh_o    = SubClass_1_1.WrapUp_msh(msh_o);
            %  > Compute PDE error.
            [X_o,Norm_o,~,~] = Class_2.Compute_ErrorPDE(msh_o);
        end
    end
end