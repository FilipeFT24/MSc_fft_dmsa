classdef B
    methods(Static)
        %% > Wrap-up B.
        function [bnd,blk] = WrapUp_B(inp,msh)
            %  > Wrap up B_1_1.
            [bnd,blk] = B_1_1.WrapUp_B_1_1(inp,msh);
        end
    end
end