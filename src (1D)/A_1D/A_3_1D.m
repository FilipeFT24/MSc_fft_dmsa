classdef A_3_1D
    methods (Static)
        %% > Wrap-up A_3 (1D).
        function [msh] = WrapUp_A_3_1D(inp,msh)
            % >> Local variables.
            np = inp.fr.np;
            
            % >> 1.
            msh = A_3_1_1D.WrapUp_A_3_1_1D(msh);
            % >> 2.
            msh = A_3_2_1D.WrapUp_A_3_2_1D(msh,np);
        end
    end
end