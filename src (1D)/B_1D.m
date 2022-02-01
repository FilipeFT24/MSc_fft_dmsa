classdef B_1D
    methods(Static)
        %% > Wrap-up B (1D).
        function [X,Norm] = WrapUp_B_1D(inp,msh)
            % >> Local variables.
            Xv_i = inp.msh.lim.Xv_i;
            Xv_f = inp.msh.lim.Xv_f;
            ft   = inp.fr.ft;
            st   = inp.fr.st;
            np   = inp.fr.np;
            ng   = inp.fr.ng;
            v    = inp.pr.v;
            g    = inp.pr.g;

            %  > Set analytic functions/values.
            pde = B_1_1_1D.WrapUp_B_1_1_1D(msh,Xv_i,Xv_f,v,g,'exp');
            
            [X,Norm] = SubClass_2_2.WrapUp_2_2(inp,msh,pde,ft,st,ng,np,v,g);
            %  > Solve...
           %pde = B_2_2_1D.WrapUp_B_2_2();
        end
    end
end