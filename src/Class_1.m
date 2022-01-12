classdef Class_1
    methods(Static)
        %% > Wrap up Class 1.
        function [inp,msh] = Set_Inputs()
            %  > Wrap up SubClass_1_1.
            [inp,msh] = SubClass_1_1.WrapUp_1_1();
            %  > Wrap up SubClass_1_2.
            msh       = SubClass_1_2.WrapUp_1_2(inp,msh); 
        end
    end
end