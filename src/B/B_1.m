classdef B_1
    methods (Static)
        %% > Wrap-up B_1.
        function [pde] = WrapUp_B_1(inp,msh)
            % >> Local variables.
            vx = inp.pr.vx;
            vy = inp.pr.vy;
            gx = inp.pr.gx;
            gy = inp.pr.gy;
            ng = inp.fr.ng;
            
            % >> 1.
            pde = B_1_1.WrapUp_B_1_1(msh,vx,vy,gx,gy);
            % >> 2.
            pde = B_1_2.WrapUp_B_1_2(msh,pde,ng);
        end            
    end
end