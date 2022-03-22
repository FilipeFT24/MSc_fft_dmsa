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
                %  > #1.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_V1_1_1D.Plot_1_1(obj,msh,fig,[exp,0]);
                subplot(1,2,2);
                Fig_V1_1_1D.Plot_1_2(obj,msh,fig,[exp,0]);
                %  > #2.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_V1_1_1D.Plot_2_1(obj,msh,fig,[exp,0]);
                subplot(1,2,2);
                Fig_V1_1_1D.Plot_2_2(obj,msh,fig,[exp,0]);
            else
                %  > #1.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V1_1_1D.Plot_1_1(obj,msh,fig,[exp,0]);
                %  > #2.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V1_1_1D.Plot_1_2(obj,msh,fig,[exp,0]);
                %  > #3.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V1_1_1D.Plot_2_1(obj,msh,fig,[exp,0]);
                %  > #4.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V1_1_1D.Plot_2_2(obj,msh,fig,[exp,0]);
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1_1(obj,msh,fig,edt)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            fig.M = repelem(":o",3);

            %  > Plot variables.
            %  > #1.
            n          = 1;
            L1{1}      = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla}\phi}|$";
            L1{2}      = "$|\bar{\tau}_{f^{\left(p\right)}}^{\nabla\phi}|$";
            L1{3}      = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla\phi}}|$";
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,L1,msh.f.Xv,obj.e.p{n}.t.f_abs);
            %  > #2.
            n          = 1;
            L2{1}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla}\phi}\|_{1}$";
            L2{2}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\nabla\phi}\|_{1}$";
            L2{3}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla\phi}}\|_{1}$";
            [L2,P2,Y2] = Fig_Tools_1D.Var_2(fig,L2,msh.f.Xv,obj.e.p{n}.t.n_abs.f);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,[L1,L2],[P1,P2],msh.f.Xv,[Y1;Y2],[-1,1],2); 
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(edt,"Fig_1_1",fig.Folder);
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_1_2(obj,msh,fig,edt)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            fig.M = repelem(":o",2);
   
            %  > Plot variables.
            %  > #1.
            n          = 1;
            L1{1}      = "$|e_{c^{\left(p\right)}}|$";
            L1{2}      = "$|\bar{\tau}_{c^{\left(p\right)}}|$";
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,L1,msh.c.Xc,[obj.e.p{n}.c.c_abs,obj.e.p{n}.t.c_abs]);
            %  > #2.
            n          = 1;
            L2{1}      = "$\|e_{c^{\left(p\right)}}\|_{1}$";
            L2{2}      = "$\|\bar{\tau}_{c^{\left(p\right)}}\|_{1}$";
            [L2,P2,Y2] = Fig_Tools_1D.Var_2(fig,L2,msh.c.Xc,[obj.e.p{n}.c.n_abs,obj.e.p{n}.t.n_abs.c]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,[L1,L2],[P1,P2],msh.f.Xv,[Y1;Y2],[-1,1],2);  
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(edt,"Fig_1_2",fig.Folder);
        end
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        function [] = Plot_2_1(obj,msh,fig,edt)
            %  > Auxiliary variables.
            Flag  = 0;
            fig.C = linspecer(9,'qualitative');

            %  > Plot variables.
            if ~Flag
                i          = 1:2;
                %  > #1.
                fig.M      = repelem(":o",2);
                L1{1}      = "$|\bar{\tau}_{f}^{\phantom{\nabla}\phi}|$";
                L1{2}      = "$|\bar{\tau}_{f}^{\nabla\phi}|$";
                [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,L1,msh.f.Xv,obj.e.a.t.f_abs(:,i));
                %  > #2.
                fig.M      = repelem("-v",2);
                n          = 1;
                L2{1}      = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\phi}}|$";
                L2{2}      = "$|\bar{\tau}_{f^{\left(p\right)}}^{\nabla\phi}|$";
                [L2,P2,Y2] = Fig_Tools_1D.Var_1(fig,L2,msh.f.Xv,obj.e.p{n}.t.f_abs(:,i));
            else
                i          = 1:3;
                %  > #1.
                fig.M      = repelem(":o",3);
                L1{1}      = "$|\bar{\tau}_{f}^{\phantom{\nabla}\phi}|$";
                L1{2}      = "$|\bar{\tau}_{f}^{\nabla\phi}|$";
                L1{3}      = "$|\bar{\tau}_{f}^{\phantom{\nabla\phi}}|$";
                [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,L1,msh.f.Xv,obj.e.a.t.f_abs(:,i)); 
                %  > #2.
                fig.M      = repelem("-v",3);
                n          = 1;
                L2{1}      = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla}\phi}|$";
                L2{2}      = "$|\bar{\tau}_{f^{\left(p\right)}}^{\nabla\phi}|$";
                L2{3}      = "$|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla\phi}}|$";
                [L2,P2,Y2] = Fig_Tools_1D.Var_1(fig,L2,msh.f.Xv,obj.e.p{n}.t.f_abs(:,i));        
            end
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,[L1,L2],[P1,P2],msh.f.Xv,[Y1;Y2],[-1,1],2); 
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(edt,"Fig_1_3",fig.Folder);
        end
        % >> 3.2. ---------------------------------------------------------
        function [] = Plot_2_2(obj,msh,fig,edt)
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C = linspecer(9,'qualitative');
            
            %  > Plot variables.
            %  > #1.
            fig.M      = repelem(":o",2);
            L1{1}      = "$|e_{c}|$";
            L1{2}      = "$|\bar{\tau}_{c}|$";
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,L1,msh.c.Xc,[obj.e.a.c.c_abs,obj.e.a.t.c_abs]);
            %  > #2.
            fig.M      = repelem("-v",2);
            n          = 1;
            L2{1}      = "$|e_{c^{\left(p\right)}}|$";
            L2{2}      = "$|\bar{\tau}_{c^{\left(p\right)}}|$";
            [L2,P2,Y2] = Fig_Tools_1D.Var_1(fig,L2,msh.c.Xc,[obj.e.p{n}.c.c_abs,obj.e.p{n}.t.c_abs]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,[L1,L2],[P1,P2],msh.f.Xv,[Y1;Y2],[-1,1],2); 
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(edt,"Fig_1_4",fig.Folder);
        end
    end
end