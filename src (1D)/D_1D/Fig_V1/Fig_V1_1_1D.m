classdef Fig_V1_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(obj,msh)
            %  > Auxiliary variables.
            exp  = 0;
            fig  = Fig_Tools_1D.Set_fig(0,exp);
            j    = 1;
            n    = 1;
            plot = [0,1];
            
            if ~exp
                %  > #1.
                if plot(1)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    Fig_V1_1_1D.Plot_1_1(obj,msh,fig,j,n,0);
                    subplot(1,2,2);
                    Fig_V1_1_1D.Plot_1_2(obj,msh,fig,j,n,0);
                end
                %  > #2.
                if plot(2)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    Fig_V1_1_1D.Plot_2_1(obj,msh,fig,j,n,0);
                    subplot(1,2,2);
                    Fig_V1_1_1D.Plot_2_2(obj,msh,fig,j,n,0);
                end
            else
                %  > #1.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V1_1_1D.Plot_1_1(obj,msh,fig,j,n,0);
                %  > #2.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V1_1_1D.Plot_1_2(obj,msh,fig,j,n,0);
                %  > #3.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V1_1_1D.Plot_2_1(obj,msh,fig,j,n,0);
                %  > #4.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V1_1_1D.Plot_2_2(obj,msh,fig,j,n,0);
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1_1(obj,msh,fig,j,n,zoom)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            
            %  > Plot variables.
            %  > #1.
            M1         = repelem("--o",3);
            L1{1}      = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla}\phi}|$";
            L1{2}      = "$|\bar{\tau}_{f^{\left(p\right)}}^{\nabla\phi}|$";
            L1{3}      = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla\phi}}|$";
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,msh.f.Xv,obj.e.p{n}.t.f_abs);
            %  > #2.
            M2         = repelem("-",3);
            L2{1}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla}\phi}\|_{1}$";
            L2{2}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\nabla\phi}\|_{1}$";
            L2{3}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla\phi}}\|_{1}$";
            [L2,P2,Y2] = Fig_Tools_1D.Var_2(fig,M2,L2,msh.f.Xv,obj.e.p{n}.t.n_abs.f(j,:));
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,[L1,L2],[P1,P2],msh.f.Xv,[Y1;Y2],[-1,4],2); 
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,"Fig_V1_1_1");
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_1_2(obj,msh,fig,j,n,zoom)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');

            %  > Plot variables.
            %  > #1.
            M1         = repelem("--o",2);
            L1{1}      = "$|e_{c^{\left(p\right)}}|$";
            L1{2}      = "$|\bar{\tau}_{c^{\left(p\right)}}|$";
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,msh.c.Xc,[obj.e.p{n}.c.c_abs,obj.e.p{n}.t.c_abs]);
            %  > #2.
            M2         = repelem("-",2);
            L2{1}      = "$\|e_{c^{\left(p\right)}}\|_{1}$";
            L2{2}      = "$\|\bar{\tau}_{c^{\left(p\right)}}\|_{1}$";
            [L2,P2,Y2] = Fig_Tools_1D.Var_2(fig,M2,L2,msh.c.Xc,[obj.e.p{n}.c.n_abs(j),obj.e.p{n}.t.n_abs.c(j)]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,[L1,L2],[P1,P2],msh.f.Xv,[Y1;Y2],[-1,4],2);  
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,"Fig_V1_1_2");
        end
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        function [] = Plot_2_1(obj,msh,fig,j,n,zoom)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            
            %  > Plot variables.
            %  > #1.
            M1         = ["--o",":v","-.d","--o",":v","-.d"];
            L1{1}      = "$|\bar{\tau}_{f^{\left(a\right)\phantom{-p}}}^{\phantom{\nabla}\phi}|$";
            L1{2}      = "$|\bar{\tau}_{f^{\left(p\right)\phantom{-p}}}^{\phantom{\nabla}\phi}|$";
            L1{3}      = "$|\bar{\tau}_{f^{\left(a-p\right)}}^{\phantom{\nabla}\phi}|$";
            L1{4}      = "$|\bar{\tau}_{f^{\left(a\right)\phantom{-p}}}^{\nabla\phi}|$";
            L1{5}      = "$|\bar{\tau}_{f^{\left(p\right)\phantom{-p}}}^{\nabla\phi}|$";
            L1{6}      = "$|\bar{\tau}_{f^{\left(a-p\right)}}^{\nabla\phi}|$";
            Xv         = repelem(msh.f.Xv,1,6);
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,Xv,...
                [obj.e.a{n}.t.f_abs(:,1),obj.e.p{n}.t.f_abs(:,1),obj.e.d{n}.t.f_abs(:,1),...
                 obj.e.a{n}.t.f_abs(:,2),obj.e.p{n}.t.f_abs(:,2),obj.e.d{n}.t.f_abs(:,2)]);
            %  > #2.
            M2         = repelem("-",6);
            L2{1}      = "$\|\bar{\tau}_{f^{\left(a\right)\phantom{-p}}}^{\phantom{\nabla}\phi}\|_{1}$";
            L2{2}      = "$\|\bar{\tau}_{f^{\left(p\right)\phantom{-p}}}^{\phantom{\nabla}\phi}\|_{1}$";
            L2{3}      = "$\|\bar{\tau}_{f^{\left(a-p\right)}}^{\phantom{\phi}}\|_{1}$";
            L2{4}      = "$\|\bar{\tau}_{f^{\left(a\right)\phantom{-p}}}^{\nabla\phi}\|_{1}$";
            L2{5}      = "$\|\bar{\tau}_{f^{\left(p\right)\phantom{-p}}}^{\nabla\phi}\|_{1}$";
            L2{6}      = "$\|\bar{\tau}_{f^{\left(a-p\right)}}^{\nabla\phi}\|_{1}$";
            Xv         = repelem(msh.f.Xv,1,6);
            [L2,P2,Y2] = Fig_Tools_1D.Var_2(fig,M2,L2,Xv,...
                [obj.e.a{n}.t.n_abs.f(j,1),obj.e.p{n}.t.n_abs.f(j,1),obj.e.d{n}.t.n_abs.f(j,1),...
                 obj.e.a{n}.t.n_abs.f(j,2),obj.e.p{n}.t.n_abs.f(j,2),obj.e.d{n}.t.n_abs.f(j,2)]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,[L1,L2],[P1,P2],msh.f.Xv,[Y1;Y2],[-1,4],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,"Fig_V1_1_3");
        end
        % >> 3.2. ---------------------------------------------------------
        function [] = Plot_2_2(obj,msh,fig,j,n,zoom)
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C = linspecer(9,'qualitative');
            
            %  > Plot variables.
            %  > #1.
            M1         = ["--o","--o",":v",":v"];
            L1{1}      = "$|e_{c^{\left(a\right)}}|$";
            L1{2}      = "$|e_{c^{\left(p\right)}}|$";
            L1{3}      = "$|\bar{\tau}_{c^{\left(a\right)}}|$";
            L1{4}      = "$|\bar{\tau}_{c^{\left(p\right)}}|$";
            Xc         = repelem(msh.c.Xc,1,4);
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,Xc,...
                [obj.e.a{n}.c.c_abs,obj.e.p{n}.c.c_abs,obj.e.a{n}.t.c_abs,obj.e.p{n}.t.c_abs]);
            %  > #2.
            M2         = repelem("-",4);
            L2{1}      = "$\|e_{c^{\left(a\right)}}\|_{1}$";
            L2{2}      = "$\|e_{c^{\left(p\right)}}\|_{1}$";
            L2{3}      = "$\|\bar{\tau}_{c^{\left(a\right)}}\|_{1}$";
            L2{4}      = "$\|\bar{\tau}_{c^{\left(p\right)}}\|_{1}$";
            Xc         = repelem(msh.c.Xc,1,4);
            [L2,P2,Y2] = Fig_Tools_1D.Var_2(fig,M2,L2,Xc,...
                [obj.e.a{n}.c.n_abs(1),obj.e.p{n}.c.n_abs(j),obj.e.a{n}.t.n_abs.c(1),obj.e.p{n}.t.n_abs.c(j)]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,[L1,L2],[P1,P2],msh.f.Xv,[Y1;Y2],[-1,4],2); 
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,"Fig_V1_1_4");
        end
    end
end