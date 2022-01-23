classdef B
    methods(Static)
        %% > Wrap-up B.
        function [pde] = WrapUp_B(inp,msh)
            % >> 1.
            pde = B_1.WrapUp_B_1(inp,msh);
            % >> 2.
            pde = B_2.WrapUp_B_2(inp,msh,pde);
        end
    end
end