classdef B
    methods(Static)
        %% > Wrap-up B.
        function [pde] = WrapUp_B(inp,msh)
            % >> Local variables.
            Xv_i = inp.msh.lim.Xv_i;
            Xv_f = inp.msh.lim.Xv_f;
            Yv_i = inp.msh.lim.Yv_i;
            Yv_f = inp.msh.lim.Yv_f;
            vx   = inp.pr.vx;
            vy   = inp.pr.vy;
            gx   = inp.pr.gx;
            gy   = inp.pr.gy;
            np_x = inp.fr.np_x;
            np_y = inp.fr.np_y;
            ng   = inp.fr.ng;
            ft   = inp.fr.ft;
            st   = inp.fr.st;
            wf   = inp.fr.wf;
            
            %  > Auxiliary arrays.
            j         = 1:size(msh.bnd.f,2);
            bnd_ff(j) = [msh.bnd.f{2,j}];
            bnd_fc(j) = [msh.bnd.f{3,j}];

            %  > Set analytic functions/values.
            pde = B_1_1.WrapUp_B_1_1(msh,Xv_i,Xv_f,Yv_i,Yv_f,vx,vy,gx,gy,'sin');
            %  > Solve...
            pde = B_2_2.WrapUp_B_2_2(msh,pde,np_x,np_y,ng,wf,vx,vy,gx,gy,bnd_fc,bnd_ff,ft,st);
        end
    end
end