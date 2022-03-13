classdef Fig_1_1D
    methods (Static)
        function [] = Plot(msh,pde)
            %  > Auxiliary variables.
            Exp = 0;
            fig = Fig_Tools_1D.Set_fig(Exp);
            N   = [1,2];

            if ~Exp
                figure(N(1)); set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_1_1D.Plot_1(msh,pde,fig,Exp,N(1));
                subplot(1,2,2);
                Fig_1_1D.Plot_2(msh,pde,fig,Exp,N(2));
            else
                figure(N(1)); set(gcf,'Units','pixels','Position',fig.Position);
                Fig_1_1D.Plot_1(msh,pde,fig,Exp,N(1));
                figure(N(2)); set(gcf,'Units','pixels','Position',fig.Position);
                Fig_1_1D.Plot_2(msh,pde,fig,Exp,N(2));
            end
        end
        % >> 1. -----------------------------------------------------------
        function [] = Plot_1(msh,pde,fig,Exp,N)
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C     = linspecer(9,'qualitative');
            fig.M     = repelem(":o",size(pde.e.t.f_abs,2));
            fig.L1{1} = "$|\bar{\tau}_{f}^{\phantom{\nabla}\phi}|$";
            fig.L1{2} = "$|\bar{\tau}_{f}^{\nabla\phi}|$";
            fig.L1{3} = "$|\bar{\tau}_{f}^{\phantom{\nabla\phi}}|$";
            fig.L1{4} = "$\|\bar{\tau}_{f}^{\phantom{\nabla}\phi}\|_{1}$";
            fig.L1{5} = "$\|\bar{\tau}_{f}^{\nabla\phi}\|_{1}$";
            fig.L1{6} = "$\|\bar{\tau}_{f}^{\phantom{\nabla\phi}}\|_{1}$";
            fig.L2{1} = "$x$";
            fig.L2{2} = "$\textrm{Error magnitude}$";
             
            %  > Plot variables.
            [~,P1,Y1] = Fig_Tools_1D.Var_1(fig,msh.f.Xv,pde.e.t.f_abs);
            [~,P2,Y2] = Fig_Tools_1D.Var_2(fig,msh.f.Xv,pde.e.t.n_abs.f);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,[P1,P2],[Y1;Y2]); 
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(N)]),fig.Folder)
            end
        end
        % >> 2. -----------------------------------------------------------
        function [] = Plot_2(msh,pde,fig,Exp,N)
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C     = linspecer(9,'qualitative');
            fig.M     = repelem(":o",2);
            fig.L1{1} = "$|e_{c}|$";
            fig.L1{2} = "$|\bar{\tau}_{c}|$";
            fig.L1{3} = "$\|e_{c}\|_{1}$";
            fig.L1{4} = "$\|\bar{\tau}_{c}\|_{1}$";
            fig.L2{1} = "$x$";
            fig.L2{2} = "$\textrm{Error magnitude}$";

            %  > Plot variables.
            [~,P1,Y1] = Fig_Tools_1D.Var_1(fig,msh.c.Xc,[pde.e.c.c_abs,pde.e.t.c_abs]);
            [~,P2,Y2] = Fig_Tools_1D.Var_2(fig,msh.c.Xc,[pde.e.c.n_abs,pde.e.t.n_abs.c]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,[P1,P2],[Y1;Y2]); 
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(N)]),fig.Folder);
            end
        end
    end
end