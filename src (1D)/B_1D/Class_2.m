classdef Class_2
    methods(Static)
        %% > Wrap up Class 2.
        function [X,Norm,bnd,blk] = Compute_ErrorPDE(inp,msh)
            %  > Wrap up SubClass_2_1.
            [fn,bnd,blk] = SubClass_2_1.WrapUp_2_1(inp,msh);
            %  > Wrap up SubClass_2_2.
            [X,Norm] = SubClass_2_2.WrapUp_2_2(inp,msh,fn,bnd,blk);
        end
    end
end