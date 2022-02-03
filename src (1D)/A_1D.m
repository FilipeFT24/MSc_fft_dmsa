classdef A_1D
    methods(Static)
        %% > Wrap-up A (1D).
        function [inp,msh] = WrapUp_A_1D(h)
            [inp]     = A_1_1D.Set_inp(h);
            [inp,msh] = A_2_1D.WrapUp_A_2_1D(inp);
        end
    end
end