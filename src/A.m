classdef A
    methods(Static)
        %% > Wrap-up A.
        % >> --------------------------------------------------------------
        % >> 1. Wrap-up A_1.
        % >> 2. Wrap-up A_2.
        % >> 3. Wrap-up A_3.
        % >> 4. Sort 'msh' fields.
        % >> --------------------------------------------------------------
        function [inp,msh] = WrapUp_A(h)
            %  > 1.
            inp = A_1.WrapUp_A_1(h);
            %  > 2.
            msh = A_2.WrapUp_A_2(inp);
            %  > 3.
            msh = A_3.WrapUp_A_3(inp);
            %  > 4.
            msh = A_Tools.Sort_msh(msh);
        end
    end
end