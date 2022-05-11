classdef Fig_2D_3
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(msh,obj)
            %  > Auxiliary variables.
            exp  = 0;
            run  = 1;
            zoom = 0;
            fig  = Fig_Tools.Set_fig(exp,run,zoom);
            x.a  = 1;
            x.b  = 1;
            
            if ~exp
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_2D_3.Plot_1(msh,obj,x,fig);
            else
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_2D_3.Plot_1(msh,obj,x,fig);
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1(msh,obj,x,fig)
            %  > Auxiliary variables.
            fig.fid = "2D_3";
            j       = x.a;
            n       = x.b;

            %  > Select variables.
            for i = 1:size(obj,2)
                h     (i,1) = msh(i).d.h;
                NNZ   (i,1) = obj(i).m  {n}.nnz.At;
                V  {1}(i,:) = obj(i).e.a{n}.t.n_abs.f(j,:);
                V  {2}(i,1) = obj(i).e.a{n}.c.n_abs  (j);
                V  {3}(i,1) = obj(i).e.a{n}.t.n_abs.c(j);
            end
            %  > Plot variables.
            L        = Fig_2D_3.Set_Legend(j);
            M        = ["-v","-^","-d","-.o","-.s"];
            [L,P1,Y] = Fig_Tools.Var_1D_1(fig,M,L,NNZ,[V{1},V{2},V{3}]);
            %  > Axis/legend,etc.
            Fig_Tools.Map_1D_2(fig,L,P1,0,NNZ,Y,[-0.5,0],2);
        end
        % >> 2.2 -------------------------------------------------------
        function [L] = Set_Legend(j)
            S(1) = Fig_Tools.Set_str_1(1);
            S(2) = Fig_Tools.Set_str_1(2);
            S(3) = Fig_Tools.Set_str_3(j);
            L{1} = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
            L{2} = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
            L{3} = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{\phantom{\nabla\phi}}\|_{_{",S(3),"}}$"]);
            L{4} = join(["$\|e_{c^{\left(a\right)}}\|_{_{",S(3),"}}$"]);
            L{5} = join(["$\|\bar{\tau}_{c^{\left(a\right)}}\|_{_{",S(3),"}}$"]);
        end
    end
end