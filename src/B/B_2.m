classdef B_2
    methods (Static)
        %% > Wrap-up B_2.
        function [pde] = WrapUp_B_2(inp,msh,pde)
            % >> Local variables.
            vx  = inp.pr.vx;
            vy  = inp.pr.vy;
            gx  = inp.pr.gx;
            gy  = inp.pr.gy;
            np  = inp.fr.np;
            ng  = inp.fr.ng;
            ft  = inp.fr.ft;
            st  = inp.fr.st;
            wf  = inp.fr.wf;
            
            %  > Auxiliary arrays.
            j         = 1:size(msh.bnd.f,2);
            bnd_ff(j) = [msh.bnd.f{2,j}];
            bnd_fc(j) = [msh.bnd.f{3,j}];
            
            % >> 1.
            pde = B_2_1.WrapUp_B_2_1(pde,np,wf);
            % >> 2.
            pde = B_2_2.WrapUp_B_2_2(msh,pde,vx,vy,gx,gy,ft,st,wf,bnd_ff,bnd_fc,ng);
        end            
    end
end