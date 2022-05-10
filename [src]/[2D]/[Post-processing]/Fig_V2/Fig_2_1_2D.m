classdef Fig_2_1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(msh,obj)
            %  > Auxiliary variables.
            exp = 0;
            fig = Fig_Tools_2D.Set_fig(1,exp);
            x.a = 1;
            x.b = 1;
            
            if ~exp
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_2_1_2D.Plot_1_1(msh,obj,x,fig,0);
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1_1(msh,obj,x,fig,zoom)
            %  > Auxiliary variables.
            fid  = "V1_0_2D";
            j    = x.a;
            n    = x.b;

            %  > Set variables to plot...
            for i = 1:size(obj,2)
                h   (i,1) = msh(i).d.h;
                NNZ (i,1) = obj(i).m  {n}.nnz.At;
                V{1}(i,1) = obj(i).e.a{n}.c.n_abs  (j);
                V{2}(i,1) = obj(i).e.a{n}.t.n_abs.c(j);
                V{3}(i,:) = obj(i).e.a{n}.t.n_abs.f(j,:);
            end
            %  > Plot variables.
            M1         = ["--o",":v","--o",":v"];
            L1{1}      = join(["$\|\bar{\tau}_{c^{\left(a\right)}}\|_{1}$"]);
            L1{2}      = join(["$\|\bar{\tau}_{c^{\left(p\right)}}\|_{1}$"]);
            L1{3}      = join(["$\|e_{c^{\left(a\right)}}\|_{1}$"]);
            L1{4}      = join(["$\|e_{c^{\left(p\right)}}\|_{1}$"]);
        end
    end
end