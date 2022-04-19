classdef Fig_V2_3_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(V)
            %  > Auxiliary variables.
            exp = 0;
            fig = Fig_Tools_1D.Set_fig(1,exp);
            j   = 1;
            n   = 1;
            
            if ~exp
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_V2_3_1D.Plot_1_1(V.inp,V.msh,V.obj,fig,j,n,0);
                subplot(1,2,2);
                Fig_V2_3_1D.Plot_1_2(V.inp,V.msh,V.obj,fig,j,n,0);
            else
                %  > #1.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V2_3_1D.Plot_1_1(V.inp,V.msh,V.obj,fig,j,n,0);
                %  > #2.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V2_3_1D.Plot_1_2(V.inp,V.msh,V.obj,fig,j,n,0);
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1_1(inp,msh,obj,fig,j,n,zoom)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            str   = Fig_Tools_1D.Set_str_1(j);
            p     = inp.ps.p+1;
                        
            %  > Set variables to plot...
            k = 1:2;
            for i = 1:size(obj,2)
                h   (i,1) = msh(i).d.h;
                NNZ (i,1) = obj(i).m.nnz.At;
                V{1}(i,k) = obj(i).e.a{n}.t.n_abs.f(j,k);
                V{2}(i,k) = obj(i).e.p{n}.t.n_abs.f(j,k);
                V{3}(i,k) = obj(i).e.d{n}.t.n_abs.f(j,k); 
            end
            for i = 1:size(V,2) 
                s{i} = C_1D.Slope(h,V{i});
            end
            %  > Plot variables.
            M1         = ["--o",":v","-.d","--o",":v","-.d"];
            L1{1}      = join(["$\|\bar{\tau}_{f^{\left(a\right)\phantom{-p}}}^{\phantom{\nabla}\phi}\|_{",str,"}$"]);
            L1{2}      = join(["$\|\bar{\tau}_{f^{\left(p\right)\phantom{-p}}}^{\phantom{\nabla}\phi}\|_{",str,"}$"]);
            L1{3}      = join(["$\|\bar{\tau}_{f^{\left(a-p\right)}}^{\phantom{\nabla}\phi}\|_{",str,"}$"]);
            L1{4}      = join(["$\|\bar{\tau}_{f^{\left(a\right)\phantom{-p}}}^{\nabla\phi}\|_{",str,"}$"]);
            L1{5}      = join(["$\|\bar{\tau}_{f^{\left(p\right)\phantom{-p}}}^{\nabla\phi}\|_{",str,"}$"]);
            L1{6}      = join(["$\|\bar{\tau}_{f^{\left(a-p\right)}}^{\nabla\phi}\|_{",str,"}$"]);
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,NNZ,...
                [V{1}(:,1),V{2}(:,1),V{3}(:,1),V{1}(:,2),V{2}(:,2),V{3}(:,2)]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,L1,P1,NNZ,Y1,[-1,1],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,"Fig_V2_1_1");
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_1_2(inp,msh,obj,fig,j,n,zoom)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            str   = Fig_Tools_1D.Set_str_1(j);
            p     = inp.ps.p+1;
            
            %  > Set variables to plot...
            for i = 1:size(obj,2)
                h   (i,1) = msh(i).d.h;
                NNZ (i,1) = obj(i).m.nnz.At;
                V{1}(i,1) = obj(i).e.a{n}.c.n_abs  (j);
                V{2}(i,1) = obj(i).e.a{n}.t.n_abs.c(j);
                V{3}(i,1) = obj(i).e.p{n}.c.n_abs  (j);
                V{4}(i,1) = obj(i).e.p{n}.t.n_abs.c(j);
            end
            for i = 1:size(V,2) 
                s(:,i) = C_1D.Slope(h,V{i});
            end
            %  > Plot variables.
            M1         = ["--o",":v","--o",":v"];
            L1{1}      = join(["$\|\bar{\tau}_{c^{\left(a\right)}}\|_{",str,"}$"]);
            L1{2}      = join(["$\|\bar{\tau}_{c^{\left(p\right)}}\|_{",str,"}$"]);
            L1{3}      = join(["$\|e_{c^{\left(a\right)}}\|_{",str,"}$"]);
            L1{4}      = join(["$\|e_{c^{\left(p\right)}}\|_{",str,"}$"]);
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,NNZ,[V{2},V{4},V{1},V{3}]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,L1,P1,NNZ,Y1,[-1,1],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,"Fig_V2_1_2");
        end
    end
end