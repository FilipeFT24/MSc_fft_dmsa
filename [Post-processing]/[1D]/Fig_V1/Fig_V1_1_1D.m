classdef Fig_V1_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,msh,obj,plot)
            %  > Auxiliary variables.
            exp = 0;
            fig = Fig_Tools_1D.Set_fig(0,exp);
            x.a = 1; %  > (:,j).
            x.b = 1; %  >   {n}.
            
            if ~exp
                if plot
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    Fig_V1_1_1D.Plot_1_1(msh.f.Xv,obj.e,x,fig,0);
                    subplot(1,2,2);
                    Fig_V1_1_1D.Plot_1_2(msh.c.Xc,msh.f.Xv,obj.e,x,fig,0);
                end
            else
                if plot
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    Fig_V1_1_1D.Plot_1_1(msh.f.Xv,obj.e,x,fig,0);
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    Fig_V1_1_1D.Plot_1_2(msh.c.Xc,msh.f.Xv,obj.e,x,fig,0);
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 2.1.1. -------------------------------------------------------
        function [] = Plot_1_1(Xv,obj_e,x,fig,zoom)
            %  > Auxiliary variables.
            fid  = "V1_1_1D (1)";
            j    = x.a;
            n    = x.b;
            plot = [1,0];
            L1   = Fig_V1_1_1D.Set_Legend_1(plot,0);
            L2   = Fig_V1_1_1D.Set_Legend_1(plot,x.a);
            
            %  > Select variables.
            if ~plot(1)
                %  > #1 (Error distribution).
                M1 = repelem(":v",2);
                V1 = obj_e.p{n}.t.f_abs  (:,1:2);
                %  > #2 (Error norms).
                M2 = repelem("-",2);
                V2 = obj_e.p{n}.t.n_abs.f(j,1:2);
            else
                if ~plot(2)
                    %  > #1 (Error distribution).
                    M1 = ["--o",":v","-.d","--o",":v","-.d"];
                    f1 = ["da","p","d"];
                    f2 = ["c","t"];
                    %  > #2 (Error norms).
                    M2 = repelem("-",6);
                else
                    %  > #1 (Error distribution).
                    M1 = ["--o","--s",":v","-.d","--o","--s",":v","-.d"];
                    f1 = ["a","da","p","d"];
                    f2 = ["c","t"];
                    %  > #2 (Error norms).
                    M2 = repelem("-",8);
                end
                l1 = length(f1);
                l2 = length(f2);
                for a = 1:l2
                    for b = 1:l1
                        V1(:,l1*(a-1)+b) = obj_e.(f1(b)){n}.t.f_abs  (:,a);
                        V2(  l1*(a-1)+b) = obj_e.(f1(b)){n}.t.n_abs.f(j,a);  
                    end
                end
            end
            %  > Plot variables.
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,Xv,V1);
            [L2,P2,Y2] = Fig_Tools_1D.Var_3(fig,M2,L2,Xv,V2);
            
            %  > Axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,[L1,L2],[P1,P2],1,Xv,[Y1;Y2],[-1,3],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,fid);
        end
        %  > 2.1.2. -------------------------------------------------------
        function [L] = Set_Legend_1(plot,n)
            S(1) = Fig_Tools_1D.Set_str_1(1);
            S(2) = Fig_Tools_1D.Set_str_1(2);
            switch plot(1)
                case 0
                    %  > w/o analytic.
                    switch n
                        case 0
                            L{1} = join(["$|\bar{\tau}_{f^{\left(p\right)}}^{",S(1),"}|$"]);
                            L{2} = join(["$|\bar{\tau}_{f^{\left(p\right)}}^{",S(2),"}|$"]);
                        case 1
                            S(3) = Fig_Tools_1D.Set_str_3(n);
                            L{1} = join(["$\|\bar{\tau}_{f^{\left(p\right)}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
                            L{2} = join(["$\|\bar{\tau}_{f^{\left(p\right)}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
                        otherwise
                            return;
                    end
                case 1
                    %  > w/  analytic.
                    switch n
                        case 0
                            switch plot(2)
                                case 0
                                    L{1} = join(["$|\bar{\tau}_{f^{\left(da\right)\phantom{-p}}}^{",S(1),"}|$"]);
                                    L{2} = join(["$|\bar{\tau}_{f^{\left(p\right)\phantom{-da}}}^{",S(1),"}|$"]);
                                    L{3} = join(["$|\bar{\tau}_{f^{\left(da-p\right)}}^{",S(1),"}|$"]);
                                    L{4} = join(["$|\bar{\tau}_{f^{\left(da\right)\phantom{-p}}}^{",S(2),"}|$"]);
                                    L{5} = join(["$|\bar{\tau}_{f^{\left(p\right)\phantom{-da}}}^{",S(2),"}|$"]);
                                    L{6} = join(["$|\bar{\tau}_{f^{\left(da-p\right)}}^{",S(2),"}|$"]);
                                case 1
                                    L{1} = join(["$|\bar{\tau}_{f^{\left(a\right)\phantom{-dp}}}^{",S(1),"}|$"]);
                                    L{2} = join(["$|\bar{\tau}_{f^{\left(da\right)\phantom{-p}}}^{",S(1),"}|$"]);
                                    L{3} = join(["$|\bar{\tau}_{f^{\left(p\right)\phantom{-da}}}^{",S(1),"}|$"]);
                                    L{4} = join(["$|\bar{\tau}_{f^{\left(da-p\right)}}^{",S(1),"}|$"]);
                                    L{5} = join(["$|\bar{\tau}_{f^{\left(a\right)\phantom{-dp}}}^{",S(2),"}|$"]);
                                    L{6} = join(["$|\bar{\tau}_{f^{\left(da\right)\phantom{-p}}}^{",S(2),"}|$"]);
                                    L{7} = join(["$|\bar{\tau}_{f^{\left(p\right)\phantom{-da}}}^{",S(2),"}|$"]);
                                    L{8} = join(["$|\bar{\tau}_{f^{\left(da-p\right)}}^{",S(2),"}|$"]);
                                otherwise
                                    return;
                            end 
                        case 1
                            S(3) = Fig_Tools_1D.Set_str_3(n);
                            switch plot(2)
                                case 0                               
                                    L{1} = join(["$\|\bar{\tau}_{f^{\left(da\right)\phantom{-p}}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
                                    L{2} = join(["$\|\bar{\tau}_{f^{\left(p\right)\phantom{-da}}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
                                    L{3} = join(["$\|\bar{\tau}_{f^{\left(da-p\right)}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
                                    L{4} = join(["$\|\bar{\tau}_{f^{\left(da\right)\phantom{-p}}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
                                    L{5} = join(["$\|\bar{\tau}_{f^{\left(p\right)\phantom{-da}}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
                                    L{6} = join(["$\|\bar{\tau}_{f^{\left(da-p\right)}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
                                case 1
                                    L{1} = join(["$\|\bar{\tau}_{f^{\left(a\right)\phantom{-dp}}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
                                    L{2} = join(["$\|\bar{\tau}_{f^{\left(da\right)\phantom{-p}}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
                                    L{3} = join(["$\|\bar{\tau}_{f^{\left(p\right)\phantom{-da}}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
                                    L{4} = join(["$\|\bar{\tau}_{f^{\left(da-p\right)}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
                                    L{5} = join(["$\|\bar{\tau}_{f^{\left(a\right)\phantom{-dp}}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
                                    L{6} = join(["$\|\bar{\tau}_{f^{\left(da\right)\phantom{-p}}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
                                    L{7} = join(["$\|\bar{\tau}_{f^{\left(p\right)\phantom{-da}}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
                                    L{8} = join(["$\|\bar{\tau}_{f^{\left(da-p\right)}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
                                otherwise
                                    return;
                            end
                        otherwise
                            return;
                    end
                otherwise
                    return;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        function [] = Plot_1_2(Xc,Xv,obj_e,x,fig,zoom)
            %  > Auxiliary variables.
            fid  = "V1_1_1D (2)";
            j    = x.a;
            n    = x.b;
            plot = [1,1];
            L1   = Fig_V1_1_1D.Set_Legend_2(plot,0);
            L2   = Fig_V1_1_1D.Set_Legend_2(plot,x.a);
            
            %  > Select variables.
            if ~plot(1)
                %  > #1 (Error distribution).
                M1      = repelem("--o",2);
                V1(:,1) = obj_e.p{n}.c.c_abs;
                V1(:,2) = obj_e.p{n}.t.c_abs;
                %  > #2 (Error norms).
                M2      = repelem("-",2);
                V2  (1) = obj_e.p{n}.c.n_abs  (j);
                V2  (2) = obj_e.p{n}.t.n_abs.c(j);
            else
                if ~plot(2)
                    %  > #1 (Error distribution).
                    M1  = ["--o",":v","-.d","--o",":v","-.d"];
                    f1  = ["da","p","d"];
                    f2  = ["c","t"];
                    %  > #2 (Error norms).
                    M2  = repelem("-",6);
                else
                    %  > #1 (Error distribution).
                    M1  = ["--o","--s",":v","-.d","--o","--s",":v","-.d"];
                    f1  = ["a","da","p","d"];
                    f2  = ["c","t"];
                    %  > #2 (Error norms).
                    M2  = repelem("-",8);
                end
                l1 = length(f1);
                l2 = length(f2);
                for a = 1:l2
                    for b = 1:l1
                        V1(:,l1*(a-1)+b) = obj_e.(f1(b)){n}.(f2(a)).c_abs;
                        switch a
                            case 1
                                V2(l1*(a-1)+b) = obj_e.(f1(b)){n}.(f2(a)).n_abs  (j);
                            case 2
                                V2(l1*(a-1)+b) = obj_e.(f1(b)){n}.(f2(a)).n_abs.c(j);
                        end  
                    end
                end
            end
            %  > Plot variables.
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,Xc,V1);
            [L2,P2,Y2] = Fig_Tools_1D.Var_3(fig,M2,L2,Xc,V2);
            
            %  > Axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,[L1,L2],[P1,P2],1,Xv,[Y1;Y2],[-1,2.5],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,fid);
        end
        %  > 2.2.2. -------------------------------------------------------
        function [L] = Set_Legend_2(plot,n)
            switch plot(1)
                case 0
                    %  > w/o analytic.
                    switch n
                        case 0
                            L{1} = "$|e_{c^{\left(p\right)}}|$";
                            L{2} = "$|\bar{\tau}_{c^{\left(p\right)}}|$";
                        case 1
                            S(1) = Fig_Tools_1D.Set_str_3(n);
                            L{1} = join(["$\|e_{c^{\left(p\right)}}\|_{_{",S(1),"}}$"]);
                            L{2} = join(["$\|\bar{\tau}_{c^{\left(p\right)}}\|_{_{",S(1),"}}$"]);
                        otherwise
                            return;
                    end
                case 1
                    %  > w/  analytic.
                    switch n
                        case 0
                            switch plot(2)
                                case 0
                                    L{1} = "$|e_{c^{\left(da\right)\phantom{-p}}}|$";
                                    L{2} = "$|e_{c^{\left(p\right)\phantom{-da}}}|$";
                                    L{3} = "$|e_{c^{\left(da-p\right)}}|$";
                                    L{4} = "$|\bar{\tau}_{c^{\left(da\right)\phantom{-p}}}|$";
                                    L{5} = "$|\bar{\tau}_{c^{\left(p\right)\phantom{-da}}}|$";
                                    L{6} = "$|\bar{\tau}_{c^{\left(da-p\right)}}|$";
                                case 1
                                    L{1} = "$|e_{c^{\left(a\right)\phantom{-dp}}}|$";
                                    L{2} = "$|e_{c^{\left(da\right)\phantom{-p}}}|$";
                                    L{3} = "$|e_{c^{\left(p\right)\phantom{-da}}}|$";
                                    L{4} = "$|e_{c^{\left(da-p\right)}}|$";
                                    L{5} = "$|\bar{\tau}_{c^{\left(a\right)\phantom{-dp}}}|$";
                                    L{6} = "$|\bar{\tau}_{c^{\left(da\right)\phantom{-p}}}|$";
                                    L{7} = "$|\bar{\tau}_{c^{\left(p\right)\phantom{-da}}}|$";
                                    L{8} = "$|\bar{\tau}_{c^{\left(da-p\right)}}|$";
                                otherwise
                                    return;
                            end 
                        otherwise
                            S(1) = Fig_Tools_1D.Set_str_3(n);
                            switch plot(2)
                                case 0
                                    L{1} = join(["$\|e_{c^{\left(da\right)\phantom{-p}}}\|_{_{",S(1),"}}$"]);
                                    L{2} = join(["$\|e_{c^{\left(p\right)\phantom{-da}}}\|_{_{",S(1),"}}$"]);
                                    L{3} = join(["$\|e_{c^{\left(da-p\right)}}\|_{_{",S(1),"}}$"]);
                                    L{4} = join(["$\|\bar{\tau}_{c^{\left(da\right)\phantom{-p}}}\|_{_{",S(1),"}}$"]);
                                    L{5} = join(["$\|\bar{\tau}_{c^{\left(p\right)\phantom{-da}}}\|_{_{",S(1),"}}$"]);
                                    L{6} = join(["$\|\bar{\tau}_{c^{\left(da-p\right)}}\|_{_{",S(1),"}}$"]);
                                case 1
                                    L{1} = join(["$\|e_{c^{\left(a\right)\phantom{-dp}}}\|_{_{",S(1),"}}$"]);
                                    L{2} = join(["$\|e_{c^{\left(da\right)\phantom{-p}}}\|_{_{",S(1),"}}$"]);
                                    L{3} = join(["$\|e_{c^{\left(p\right)\phantom{-da}}}\|_{_{",S(1),"}}$"]);
                                    L{4} = join(["$\|e_{c^{\left(da-p\right)}}\|_{_{",S(1),"}}$"]);
                                    L{5} = join(["$\|\bar{\tau}_{c^{\left(a\right)\phantom{-dp}}}\|_{_{",S(1),"}}$"]);
                                    L{6} = join(["$\|\bar{\tau}_{c^{\left(da\right)\phantom{-p}}}\|_{_{",S(1),"}}$"]);
                                    L{7} = join(["$\|\bar{\tau}_{c^{\left(p\right)\phantom{-da}}}\|_{_{",S(1),"}}$"]);
                                    L{8} = join(["$\|\bar{\tau}_{c^{\left(da-p\right)}}\|_{_{",S(1),"}}$"]);
                                otherwise
                                    return;
                            end
                    end
                otherwise
                    return;
            end
        end
    end
end