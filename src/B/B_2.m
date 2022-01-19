classdef B_2
    methods (Static)
        %% > Wrap-up B_2.
        function [pde] = WrapUp_B_2(inp,msh,pde)
            % >> Local variables.
            vx = inp.pr.vx;
            vy = inp.pr.vy;
            gx = inp.pr.gx;
            gy = inp.pr.gy;
            np = inp.fr.np;
            ng = inp.fr.ng;
            wf = inp.fr.wf;
            
            % >> 1.
            pde = B_2_1.WrapUp_B_2_1(pde,np,wf);
            % >> 2.
            pde = B_2_2.WrapUp_B_2_2(msh,pde,ng,vx,vy,gx,gy);
        end            
    end
end