classdef A_1D
    methods(Static)
        %% > Wrap-up A (1D).
        % >> --------------------------------------------------------------
        % >> 1. Wrap-up A_1.
        % >> 2. Wrap-up A_2.
        % >> 3. Wrap-up A_3.
        % >> 4. Sort 'msh' fields.
        % >> --------------------------------------------------------------
        function [inp,msh] = WrapUp_A_1D(h)
            [inp]     = A_1_1D.WrapUp_A_1_1D(h);
            [inp,msh] = A_2_1D.WrapUp_A_2_1D(inp);
        end
    end
end