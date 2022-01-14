classdef Class_1
    methods(Static)
        %% > Wrap-up Class 1.
        % >> --------------------------------------------------------------
        % >> 1. Wrap-up SubClass_1_1.
        % >> 2. Wrap-up SubClass_1_2.
        % >> 3. Wrap-up SubClass_1_3.
        % >> --------------------------------------------------------------
        function [inp,msh] = WrapUp_1()
            %  > Wrap-up SubClass_1_1.
            inp = SubClass_1_1.WrapUp_1_1();
            %  > Wrap-up SubClass_1_2.
            msh = SubClass_1_2.WrapUp_1_2(inp);
            %  > Wrap-up SubClass_1_3.
            msh = SubClass_1_3.WrapUp_1_3(inp,msh);
        end
    end
end