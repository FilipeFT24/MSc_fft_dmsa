classdef B_1D
    methods(Static)
        %% > Wrap-up B (1D).
        function [msh,pde] = WrapUp_B_1D(inp,msh)
            % >> Local variables.
            v       = inp.pr.v;
            g       = inp.pr.g;
            bnd_w   = inp.pr.w;
            bnd_e   = inp.pr.e;
            bnd     = [string(bnd_w),string(bnd_e)];
            p       = inp.fr.p;
            p_adapt = inp.fr.p_adapt;
            
            %  > Set analytic functions/values.
            [pde]     = B_1_1D.WrapUp_B_1_1D(msh,v,g,p,"exp");
            %  > Solve PDE.
            [msh,pde] = B_2_1D.SetUp_Problem(msh,pde,v,g,bnd,p,p_adapt);
            %  > Organize structures.
            [msh]     = Tools_1D.Sort_struct(msh);
        end
    end
end