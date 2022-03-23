classdef Fig_V1_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(obj,msh)
            %  > Auxiliary variables.
            run = 0;
            exp = 0;
            fig = Fig_Tools_1D.Set_fig(run,exp);
            N   = [1,2,3,4];
            
            if ~exp
                j    = 1;
                n    = 2;
                plot = [0,1];
                %  > #1.
                if plot(1)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    Fig_V1_1_1D.Plot_1_1(obj,msh,fig,j,n,[exp,0]);
                    subplot(1,2,2);
                    Fig_V1_1_1D.Plot_1_2(obj,msh,fig,j,n,[exp,0]);
                end
                %  > #2.
                if plot(2)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    Fig_V1_1_1D.Plot_2_1(obj,msh,fig,j,n,[exp,0]);
                    subplot(1,2,2);
                    Fig_V1_1_1D.Plot_2_2(obj,msh,fig,j,n,[exp,0]);
                end
            else
                j = 1;
                n = 1;
                %  > #1.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V1_1_1D.Plot_1_1(obj,msh,fig,j,n,[exp,0]);
                %  > #2.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V1_1_1D.Plot_1_2(obj,msh,fig,j,n,[exp,0]);
                %  > #3.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V1_1_1D.Plot_2_1(obj,msh,fig,j,n,[exp,0]);
                %  > #4.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V1_1_1D.Plot_2_2(obj,msh,fig,j,n,[exp,0]);
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1_1(obj,msh,fig,j,n,edt)
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
            Fig_Tools_1D.Set_Plot(fig,[L1,L2],[P1,P2],msh.f.Xv,[Y1;Y2],[-1,1],2); 
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(edt,"Fig_V1_1_1",fig.Folder);
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_1_2(obj,msh,fig,j,n,edt)
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
            Fig_Tools_1D.Set_Plot(fig,[L1,L2],[P1,P2],msh.f.Xv,[Y1;Y2],[-1,1],2);  
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(edt,"Fig_V1_1_2",fig.Folder);
        end
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        function [] = Plot_2_1(obj,msh,fig,j,n,edt)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            
            %  > Plot variables.
            %  > #1.
            M1         = ["--o","--o",":v",":v"];
            L1{1}      = "$|\bar{\tau}_{f}^{\phantom{\nabla}\phi}|$";
            L1{2}      = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\phi}}|$";
            L1{3}      = "$|\bar{\tau}_{f}^{\nabla\phi}|$";
            L1{4}      = "$|\bar{\tau}_{f^{\left(p\right)}}^{\nabla\phi}|$";
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,msh.f.Xv,...
                [obj.e.a{n}.t.f_abs(:,1),obj.e.p{n}.t.f_abs(:,1),obj.e.a{n}.t.f_abs(:,2),obj.e.p{n}.t.f_abs(:,2)]);
            %  > #2.
            M2         = repelem("-",4);
            L2{1}      = "$\|\bar{\tau}_{f}^{\phantom{\nabla}\phi}\|_{1}$";
            L2{2}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\phi}}\|_{1}$";
            L2{3}      = "$\|\bar{\tau}_{f}^{\nabla\phi}\|_{1}$";
            L2{4}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\nabla\phi}\|_{1}$";
            [L2,P2,Y2] = Fig_Tools_1D.Var_2(fig,M2,L2,msh.f.Xv,...
                [obj.e.a{n}.t.n_abs.f(j,1),obj.e.p{n}.t.n_abs.f(j,1),obj.e.a{n}.t.n_abs.f(j,2),obj.e.p{n}.t.n_abs.f(j,2)]);           
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,[L1,L2],[P1,P2],msh.f.Xv,[Y1;Y2],[-1,1],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(edt,"Fig_V1_1_3",fig.Folder);
        end
        % >> 3.2. ---------------------------------------------------------
        function [] = Plot_2_2(obj,msh,fig,j,n,edt)
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C = linspecer(9,'qualitative');
            
            %  > Plot variables.
            %  > #1.
            M1         = ["--o","--o",":v",":v"];
            L1{1}      = "$|e_{c\phantom{^{\left(p\right)}}}|$";
            L1{2}      = "$|e_{c^{\left(p\right)}}|$";
            L1{3}      = "$|\bar{\tau}_{c\phantom{^{\left(p\right)}}}|$";
            L1{4}      = "$|\bar{\tau}_{c^{\left(p\right)}}|$";
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,msh.c.Xc,...
                [obj.e.a{n}.c.c_abs,obj.e.p{n}.c.c_abs,obj.e.a{n}.t.c_abs,obj.e.p{n}.t.c_abs]);
            %  > #2.
            M2         = repelem("-",4);
            L2{1}      = "$\|e_{c\phantom{^{\left(p\right)}}}\|_{1}$";
            L2{2}      = "$\|e_{c^{\left(p\right)}}\|_{1}$";
            L2{3}      = "$\|\bar{\tau}_{c\phantom{^{\left(p\right)}}}\|_{1}$";
            L2{4}      = "$\|\bar{\tau}_{c^{\left(p\right)}}\|_{1}$";
            [L2,P2,Y2] = Fig_Tools_1D.Var_2(fig,M2,L2,msh.c.Xc,...
                [obj.e.a{n}.c.n_abs(1),obj.e.p{n}.c.n_abs(j),obj.e.a{n}.t.n_abs.c(1),obj.e.p{n}.t.n_abs.c(j)]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,[L1,L2],[P1,P2],msh.f.Xv,[Y1;Y2],[-1,1],2); 
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(edt,"Fig_V1_1_4",fig.Folder);
        end
    end
end