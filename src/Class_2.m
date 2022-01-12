classdef Class_2
    methods(Static)
        %% > Wrap up Class 2.
        function [msh,bnd,blk] = Compute_ErrorPDE(inp,msh)
            %  > Wrap up SubClass_2_1.
            [fn,bnd,blk] = SubClass_2_1.WrapUp_2_1(inp,msh);
            %  > Wrap up SubClass_2_2.
            msh = SubClass_2_2.WrapUp_2_2(inp,msh,fn);
        end
    end
end