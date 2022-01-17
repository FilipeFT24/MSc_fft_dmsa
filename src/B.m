classdef B
    methods(Static)
        %% > Wrap-up B.
        function [pde] = WrapUp_B(inp,msh)
            %  > Working directory.
            Tools.Set_Directory('B');
            % >> 1.
            pde = B_1_1.WrapUp_B_1_1(inp,msh);
            % >> 2.
            pde = B_2.WrapUp_B_2(inp,msh,pde);
        end
    end
end