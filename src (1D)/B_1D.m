classdef B_1D
    methods(Static)
        %% > Wrap-up B (1D).
        function [msh,pde] = WrapUp_B_1D(inp,msh)
            % >> Local variables.
            Xv_i  = inp.msh.lim.Xv_i;
            Xv_f  = inp.msh.lim.Xv_f;
            st    = inp.fr.st;
            np    = inp.fr.np;
            ng    = inp.fr.ng;
            v     = inp.pr.v;
            g     = inp.pr.g;
            bnd_w = inp.pr.w;
            bnd_e = inp.pr.e;
            bnd   = [string(bnd_w),string(bnd_e)];
            
            
            
            

            %  > Set analytic functions/values.
            [pde]     = B_1_1_1D.WrapUp_B_1_1_1D(msh,Xv_i,Xv_f,v,g,'exp');
            %  > Solve...
            [msh,pde] = B_2_1D.WrapUp_B_2_1D(msh,pde,st,ng,np,v,g,bnd);
        end
    end
end