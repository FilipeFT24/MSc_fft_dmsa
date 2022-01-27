classdef B_1
    methods (Static)
        %% > Wrap-up B_1.
        function [pde] = WrapUp_B_1(inp,msh)
            % >> Local variables.
            Xv_i = inp.msh.lim.Xv_i;
            Xv_f = inp.msh.lim.Xv_f;
            Yv_i = inp.msh.lim.Yv_i;
            Yv_f = inp.msh.lim.Yv_f;
            vx   = inp.pr.vx;
            vy   = inp.pr.vy;
            gx   = inp.pr.gx;
            gy   = inp.pr.gy;
            ng   = inp.fr.ng;
            
            % >> 1.
            pde = B_1_1.WrapUp_B_1_1(msh,Xv_i,Xv_f,Yv_i,Yv_f,vx,vy,gx,gy);
            % >> 2.
            pde = B_1_2.WrapUp_B_1_2(pde,ng);
        end            
    end
end