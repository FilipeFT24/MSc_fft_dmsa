classdef A_3_2
    methods (Static)
        %% > Wrap-up A_3_2.
        % >> --------------------------------------------------------------  
        % >> 1. Stencil setup.
        % >> 2. Stencil extension(s).
        % >> --------------------------------------------------------------
        function [msh] = WrapUp_A_3_2(msh,Type,NLay,p,et)
            %  > Auxiliary arrays.
            for i = 1:size(msh.bnd.c,2)
                bnd_cc(i) = msh.bnd.c{2,i};
            end
            
            % >> 1.
            msh = A_3_2_1.Set_Stencil(msh,bnd_cc,Type,NLay);           
            % >> 2.
            msh = A_3_2_2.Extend_Stencil(msh,bnd_cc,p,et);
        end
    end
end