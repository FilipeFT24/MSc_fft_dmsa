classdef Fig_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(msh,pde)
            %  > Auxiliary variables.
            Exp = 0;
            fig = Fig_Tools_1D.Set_fig(Exp);
            N   = [1,2,3,4];
            
            if ~Exp
                %  > #1.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_1_1D.Plot_1_1(msh,pde,fig,Exp,N(1));
                subplot(1,2,2);
                Fig_1_1D.Plot_1_2(msh,pde,fig,Exp,N(2));
                %  > #2.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_1_1D.Plot_2_1(msh,pde,fig,Exp,N(3));
                subplot(1,2,2);
                Fig_1_1D.Plot_2_2(msh,pde,fig,Exp,N(4));
            else
                %  > #1.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_1_1D.Plot_1_1(msh,pde,fig,Exp,N(1));
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_1_1D.Plot_1_2(msh,pde,fig,Exp,N(2));
                %  > #2.
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1_1(msh,pde,fig,Exp,N)
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C     = linspecer(9,'qualitative');
            fig.M     = repelem(":o",size(pde.e.t.a.f_abs,2)+size(pde.e.t.a.n_abs.f,2));
            fig.L1{1} = "$|\bar{\tau}_{f}^{\phantom{\nabla}\phi}|$";
            fig.L1{2} = "$|\bar{\tau}_{f}^{\nabla\phi}|$";
            fig.L1{3} = "$|\bar{\tau}_{f}^{\phantom{\nabla\phi}}|$";
            fig.L1{4} = "$\|\bar{\tau}_{f}^{\phantom{\nabla}\phi}\|_{1}$";
            fig.L1{5} = "$\|\bar{\tau}_{f}^{\nabla\phi}\|_{1}$";
            fig.L1{6} = "$\|\bar{\tau}_{f}^{\phantom{\nabla\phi}}\|_{1}$";
            fig.L2{1} = "$x$";
            fig.L2{2} = "$\textrm{Error magnitude}$";
             
            %  > Plot variables.
            [~,P1,Y1] = Fig_Tools_1D.Var_1(fig,msh.f.Xv,pde.e.t.a.f_abs);
            [~,P2,Y2] = Fig_Tools_1D.Var_2(fig,msh.f.Xv,pde.e.t.a.n_abs.f);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,[P1,P2],[Y1;Y2],2); 
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(N)]),fig.Folder)
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_1_2(msh,pde,fig,Exp,N)
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
            [~,P1,Y1] = Fig_Tools_1D.Var_1(fig,msh.c.Xc,[pde.e.c.c_abs,pde.e.t.a.c_abs]);
            [~,P2,Y2] = Fig_Tools_1D.Var_2(fig,msh.c.Xc,[pde.e.c.n_abs,pde.e.t.a.n_abs.c]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,[P1,P2],[Y1;Y2],2); 
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(N)]),fig.Folder);
            end
        end
        %% > 3. -----------------------------------------------------------
        function [] = Plot_2_1(msh,pde,fig,Exp,N)
            % >> 2.1. ---------------------------------------------------------
            %  > Auxiliary variables (colors/labels/etc.).
            n         = 1;
            fig.C     = linspecer(9,'qualitative');
            fig.L1{1} = "$|\bar{\tau}_{f}^{\phantom{\nabla}\phi}|$";
            fig.L1{2} = "$|\bar{\tau}_{f}^{\nabla\phi}|$";
            fig.L1{3} = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phi}|$";
            fig.L1{4} = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla}\phi}|$";
            fig.L2{1} = "$x$";
            fig.L2{2} = "$\textrm{Error magnitude}$";
             
            %  > Plot variables.
            fig.M     = repelem(":o",2);
            [~,P1,Y1] = Fig_Tools_1D.Var_1(fig,msh.f.Xv,pde.e.t.p{n}.f_abs(:,1:2));
            fig.M     = repelem("-v",2);
            [~,P2,Y2] = Fig_Tools_1D.Var_1(fig,msh.f.Xv,pde.e.t.a.f_abs   (:,1:2));
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,[P1,P2],[Y1;Y2],2); 
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(N)]),fig.Folder)
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_2_2(msh,pde,fig,Exp,N)
            %  > Auxiliary variables (colors/labels/etc.).
            n         = 1;
            fig.C     = linspecer(9,'qualitative');
            fig.L1{1} = "$|\bar{\tau}_{c}|$";
            fig.L1{2} = "$|\bar{\tau}_{c^{\left(p\right)}}|$";
            fig.L2{1} = "$x$";
            fig.L2{2} = "$\textrm{Error magnitude}$";

            %  > Plot variables.
            fig.M (1) = ":o";
            [~,P1,Y1] = Fig_Tools_1D.Var_1(fig,msh.c.Xc,pde.e.t.p{n}.c_abs);
            fig.M (1) = "-v";
            [~,P2,Y2] = Fig_Tools_1D.Var_1(fig,msh.c.Xc,pde.e.t.a.c_abs);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,[P1,P2],[Y1;Y2],2); 
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(N)]),fig.Folder);
            end
        end
    end
end