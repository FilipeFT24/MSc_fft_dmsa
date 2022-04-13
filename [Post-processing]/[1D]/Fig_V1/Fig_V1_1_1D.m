classdef Fig_V1_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,msh,obj)
            %  > Auxiliary variables.
            exp = 0;
            fig = Fig_Tools_1D.Set_fig(0,exp);
            x.a = 1; %  > (:,j).
            x.b = 1; %  >   {n}.
            
            if ~exp
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_V1_1_1D.Plot_1_1(msh.f.Xv,obj.e,x,fig,0);
                subplot(1,2,2);
                Fig_V1_1_1D.Plot_1_2(msh.c.Xc,msh.f.Xv,obj.e,x,fig,0); 
            else
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 2.1.1. -------------------------------------------------------
        function [] = Plot_1_1(Xv,obj_e,x,fig,zoom)
            %  > Auxiliary variables.
            fid   = "V1_1_1D (1)";
            j     = x.a;
            n     = x.b;
            p_all = 1;
            L1    = Fig_V1_1_1D.Set_Legend_1(p_all,0);
            L2    = Fig_V1_1_1D.Set_Legend_1(p_all,x.a);
            
            %  > Select variables.
            if ~p_all
                %  > #1 (Error distribution).
                M1 = repelem(":v",2);
                V1 = obj_e.p{n}.t.f_abs  (:,1:2);
                %  > #2 (Error norms).
                M2 = repelem("-",6);
                V2 = obj_e.p{n}.t.n_abs.f(j,1:2);
            else
                %  > #1 (Error distribution).
                M1      = ["--o",":v","-.d","--o",":v","-.d"];
                V1(:,1) = obj_e.a{n}.t.f_abs  (:,1);
                V1(:,2) = obj_e.p{n}.t.f_abs  (:,1);
                V1(:,3) = obj_e.d{n}.t.f_abs  (:,1);
                V1(:,4) = obj_e.a{n}.t.f_abs  (:,2);
                V1(:,5) = obj_e.p{n}.t.f_abs  (:,2);
                V1(:,6) = obj_e.d{n}.t.f_abs  (:,2);
                %  > #2 (Error norms).
                M2      = repelem("-",6);
                V2  (1) = obj_e.a{n}.t.n_abs.f(j,1);
                V2  (2) = obj_e.p{n}.t.n_abs.f(j,1);
                V2  (3) = obj_e.d{n}.t.n_abs.f(j,1);
                V2  (4) = obj_e.a{n}.t.n_abs.f(j,2);
                V2  (5) = obj_e.p{n}.t.n_abs.f(j,2);
                V2  (6) = obj_e.d{n}.t.n_abs.f(j,2);
            end
            %  > Plot variables.
            [L1,P1,Y1]  = Fig_Tools_1D.Var_1(fig,M1,L1,Xv,V1);
            [L2,P2,Y2]  = Fig_Tools_1D.Var_3(fig,M2,L2,Xv,V2);

            %  > Axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,[L1,L2],[P1,P2],1,Xv,[Y1;Y2],[-1,1],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,fid);
        end
        %  > 2.1.2. -------------------------------------------------------
        function [L] = Set_Legend_1(f,n)
            S(1) = Fig_Tools_1D.Set_str_1(1);
            S(2) = Fig_Tools_1D.Set_str_1(2);
            switch f
                case 0
                    %  > w/o analytic.
                    switch n
                        case 0
                            L{1} = join(["$|\bar{\tau}_{f^{\left(p\right)}}^{",S(1),"}|$"]);
                            L{2} = join(["$|\bar{\tau}_{f^{\left(p\right)}}^{",S(2),"}|$"]);
                        otherwise
                            S(3) = Fig_Tools_1D.Set_str_3(n);
                            L{1} = join(["$\|\bar{\tau}_{f^{\left(p\right)}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
                            L{2} = join(["$\|\bar{\tau}_{f^{\left(p\right)}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
                    end
                case 1
                    %  > w/  analytic.
                    switch n
                        case 0
                            L{1} = join(["$|\bar{\tau}_{f^{\left(a\right)\phantom{-p}}}^{",S(1),"}|$"]);
                            L{2} = join(["$|\bar{\tau}_{f^{\left(p\right)\phantom{-p}}}^{",S(1),"}|$"]);
                            L{3} = join(["$|\bar{\tau}_{f^{\left(a-p\right)}}^{",S(1),"}|$"]);
                            L{4} = join(["$|\bar{\tau}_{f^{\left(a\right)\phantom{-p}}}^{",S(2),"}|$"]);
                            L{5} = join(["$|\bar{\tau}_{f^{\left(p\right)\phantom{-p}}}^{",S(2),"}|$"]);
                            L{6} = join(["$|\bar{\tau}_{f^{\left(a-p\right)}}^{",S(2),"}|$"]);
                        otherwise
                            S(3) = Fig_Tools_1D.Set_str_3(n);
                            L{1} = join(["$\|\bar{\tau}_{f^{\left(a\right)\phantom{-p}}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
                            L{2} = join(["$\|\bar{\tau}_{f^{\left(p\right)\phantom{-p}}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
                            L{3} = join(["$\|\bar{\tau}_{f^{\left(a-p\right)}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
                            L{4} = join(["$\|\bar{\tau}_{f^{\left(a\right)\phantom{-p}}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
                            L{5} = join(["$\|\bar{\tau}_{f^{\left(p\right)\phantom{-p}}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
                            L{6} = join(["$\|\bar{\tau}_{f^{\left(a-p\right)}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
                    end
                otherwise
                    return;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        function [] = Plot_1_2(Xc,Xv,obj_e,x,fig,zoom)
            %  > Auxiliary variables.
            fid   = "V1_1_1D (2)";
            j     = x.a;
            n     = x.b;
            p_all = 1;
            L1    = Fig_V1_1_1D.Set_Legend_2(p_all,0);
            L2    = Fig_V1_1_1D.Set_Legend_2(p_all,x.a);
            
            %  > Select variables.
            if ~p_all
                %  > #1 (Error distribution).
                M1 = repelem("--o",2);
                V1(:,1) = obj_e.p{n}.c.c_abs;
                V1(:,2) = obj_e.p{n}.t.c_abs;
                %  > #2 (Error norms).
                M2 = repelem("-",2);
                V2  (1) = obj_e.p{n}.c.n_abs  (j);
                V2  (2) = obj_e.p{n}.t.n_abs.c(j);
            else
                %  > #1 (Error distribution).
                M1      = ["--o","--o",":v",":v"];
                V1(:,1) = obj_e.a{n}.c.c_abs;
                V1(:,2) = obj_e.p{n}.c.c_abs;
                V1(:,3) = obj_e.a{n}.t.c_abs;
                V1(:,4) = obj_e.p{n}.t.c_abs;
                %  > #2 (Error norms).
                M2      = repelem("-",4);
                V2  (1) = obj_e.a{n}.c.n_abs  (j);
                V2  (2) = obj_e.p{n}.c.n_abs  (j);
                V2  (3) = obj_e.a{n}.t.n_abs.c(j);
                V2  (4) = obj_e.p{n}.t.n_abs.c(j);
            end
            %  > Plot variables.
            [L1,P1,Y1]  = Fig_Tools_1D.Var_1(fig,M1,L1,Xc,V1);
            [L2,P2,Y2]  = Fig_Tools_1D.Var_3(fig,M2,L2,Xc,V2);

            %  > Axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,[L1,L2],[P1,P2],1,Xv,[Y1;Y2],[-1,1],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,fid);
        end
        %  > 2.2.2. -------------------------------------------------------
        function [L] = Set_Legend_2(f,n)
            switch f
                case 0
                    %  > w/o analytic.
                    switch n
                        case 0
                            L{1} = "$|e_{c^{\left(p\right)}}|$";
                            L{2} = "$|\bar{\tau}_{c^{\left(p\right)}}|$";
                        otherwise
                            S(1) = Fig_Tools_1D.Set_str_3(n);
                            L{1} = join(["$\|e_{c^{\left(p\right)}}\|_{_{",S(1),"}}$"]);
                            L{2} = join(["$\|\bar{\tau}_{c^{\left(p\right)}}\|_{_{",S(1),"}}$"]);
                    end
                case 1
                    %  > w/  analytic.
                    switch n
                        case 0
                            L{1} = "$|e_{c^{\left(a\right)}}|$";
                            L{2} = "$|e_{c^{\left(p\right)}}|$";
                            L{3} = "$|\bar{\tau}_{c^{\left(a\right)}}|$";
                            L{4} = "$|\bar{\tau}_{c^{\left(p\right)}}|$";
                        otherwise
                            S(1) = Fig_Tools_1D.Set_str_3(n);
                            L{1} = join(["$\|e_{c^{\left(a\right)}}\|_{_{",S(1),"}}$"]);
                            L{2} = join(["$\|e_{c^{\left(p\right)}}\|_{_{",S(1),"}}$"]);
                            L{3} = join(["$\|\bar{\tau}_{c^{\left(a\right)}}\|_{_{",S(1),"}}$"]);
                            L{4} = join(["$\|\bar{\tau}_{c^{\left(p\right)}}\|_{_{",S(1),"}}$"]);
                    end
                otherwise
                    return;
            end
        end   
    end
end