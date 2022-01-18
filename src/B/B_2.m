classdef B_2
    methods (Static)
        %% > Wrap-up B_2.
        function [pde] = WrapUp_B_2(inp,msh,pde)
            % >> 1.
            pde = B_2_1.WrapUp_B_2_1(inp,pde);
            % >> 2.
            pde = B_2_2.WrapUp_B_2_2(inp,msh,pde);
        end            
    end
end