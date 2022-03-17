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
                Fig_1_1D.Plot_1_1(msh,pde,fig,Exp,N(1),1);
                subplot(1,2,2);
                Fig_1_1D.Plot_1_2(msh,pde,fig,Exp,N(2),1);
                %  > #2.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_1_1D.Plot_2_1(msh,pde,fig,Exp,N(3),1);
                subplot(1,2,2);
                Fig_1_1D.Plot_2_2(msh,pde,fig,Exp,N(4),1);
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
        function [] = Plot_1_1(msh,pde,fig,Exp,N,n)
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C     = linspecer(9,'qualitative');
            fig.M     = repelem(":o",size(pde.e.a.t.f_abs,2)+size(pde.e.a.t.n_abs.f,2));
            fig.L1{1} = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla}\phi}|$";
            fig.L1{2} = "$|\bar{\tau}_{f^{\left(p\right)}}^{\nabla\phi}|$";
            fig.L1{3} = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla\phi}}|$";
            fig.L1{4} = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla}\phi}\|_{1}$";
            fig.L1{5} = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\nabla\phi}\|_{1}$";
            fig.L1{6} = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla\phi}}\|_{1}$";
            fig.L2{1} = "$x$";
            fig.L2{2} = "$\textrm{Error magnitude}$";
             
            %  > Plot variables.
            [~,P1,Y1] = Fig_Tools_1D.Var_1(fig,msh.f.Xv,pde.e.p.t{n}.f_abs);
            [~,P2,Y2] = Fig_Tools_1D.Var_2(fig,msh.f.Xv,pde.e.p.t{n}.n_abs.f);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,[P1,P2],[Y1;Y2],2); 
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(N)]),fig.Folder)
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_1_2(msh,pde,fig,Exp,N,n)
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C     = linspecer(9,'qualitative');
            fig.M     = repelem(":o",2);
            fig.L1{1} = "$|e_{c^{\left(p\right)}}|$";
            fig.L1{2} = "$|\bar{\tau}_{c^{\left(p\right)}}|$";
            fig.L1{3} = "$\|e_{c^{\left(p\right)}}\|_{1}$";
            fig.L1{4} = "$\|\bar{\tau}_{c^{\left(p\right)}}\|_{1}$";
            fig.L2{1} = "$x$";
            fig.L2{2} = "$\textrm{Error magnitude}$";

            %  > Plot variables.
            [~,P1,Y1] = Fig_Tools_1D.Var_1(fig,msh.c.Xc,[pde.e.p.c{n}.c_abs  ,pde.e.p.t{n}.c_abs]);
            [~,P2,Y2] = Fig_Tools_1D.Var_2(fig,msh.c.Xc,[pde.e.p.c{n}.n_abs.c,pde.e.p.t{n}.n_abs.c]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,[P1,P2],[Y1;Y2],2); 
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(N)]),fig.Folder);
            end
        end
        %% > 3. -----------------------------------------------------------
        function [] = Plot_2_1(msh,pde,fig,Exp,N,n)
            % >> 2.1. ---------------------------------------------------------
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C     = linspecer(9,'qualitative');
            fig.L1{1} = "$|\bar{\tau}_{f}^{\phantom{\nabla}\phi}|$";
            fig.L1{2} = "$|\bar{\tau}_{f}^{\nabla\phi}|$";
            fig.L1{3} = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phi}|$";
            fig.L1{4} = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla}\phi}|$";
            fig.L2{1} = "$x$";
            fig.L2{2} = "$\textrm{Error magnitude}$";
             
            %  > Plot variables.
            fig.M     = repelem(":o",2);
            [~,P1,Y1] = Fig_Tools_1D.Var_1(fig,msh.f.Xv,pde.e.p.t{n}.f_abs(:,1:2));
            fig.M     = repelem("-v",2);
            [~,P2,Y2] = Fig_Tools_1D.Var_1(fig,msh.f.Xv,pde.e.a.t.f_abs   (:,1:2));
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,[P1,P2],[Y1;Y2],2); 
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(N)]),fig.Folder)
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_2_2(msh,pde,fig,Exp,N,n)
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C     = linspecer(9,'qualitative');
            fig.L1{1} = "$|e_{c}|$";
            fig.L1{2} = "$|\bar{\tau}_{c}|$";
            fig.L1{3} = "$|e_{c^{\left(p\right)}}|$";
            fig.L1{4} = "$|\bar{\tau}_{c^{\left(p\right)}}|$";
            fig.L2{1} = "$x$";
            fig.L2{2} = "$\textrm{Error magnitude}$";

            %  > Plot variables.
            fig.M     = repelem(":o",2);
            [~,P1,Y1] = Fig_Tools_1D.Var_1(fig,msh.c.Xc,[pde.e.a.c.c_abs,pde.e.a.t.c_abs]);
            fig.M     = repelem("-v",2);
            [~,P2,Y2] = Fig_Tools_1D.Var_1(fig,msh.c.Xc,[pde.e.p.c{n}.c_abs,pde.e.p.t{n}.c_abs]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,[P1,P2],[Y1;Y2],2); 
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(N)]),fig.Folder);
            end
        end
    end
end