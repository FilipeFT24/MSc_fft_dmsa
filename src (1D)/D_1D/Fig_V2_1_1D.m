classdef Fig_V2_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(V)
            %  > Auxiliary variables.
            run = 1;
            exp = 0;
            fig = Fig_Tools_1D.Set_fig(run,exp);
            n   = 2;
            N   = [1,2];
            
            if ~exp
                %  > #1.
                j = 1;
                n = 2;
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_V2_1_1D.Plot_1_1(V.inp,V.obj,fig,j,n,[exp,0]);
                subplot(1,2,2);
                Fig_V2_1_1D.Plot_1_2(V.inp,V.obj,fig,j,n,[exp,0]);
            else
                %  > #1.
                j = 1;
                n = 2;
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V2_1_1D.Plot_1_1(V.inp,V.obj,fig,j,n,[exp,0]);
                %  > #2.
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_V2_1_1D.Plot_1_2(V.inp,V.obj,fig,j,n,[exp,0]);
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1_1(inp,obj,fig,j,n,edt)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            str   = Fig_Tools_1D.Set_str_1(j);
            p     = inp.ps.p+1;
                        
            %  > Set variables to plot...
            for i = 1:size(obj,2)
                NNZ(i,1) = obj(i).m.nnz.At;
                V1 (i,:) = obj(i).e.a{n}.t.n_abs.f(j,:);
                V2 (i,:) = obj(i).e.p{n}.t.n_abs.f(j,:);
            end
            %  > Plot variables.
            %  > #1.
            M1         = repelem("--o",3);
            L1{1}      = join(["$\|\bar{\tau}_{f}^{\phantom{\nabla}\phi}\|_{",str,"}$"]);
            L1{2}      = join(["$\|\bar{\tau}_{f}^{\nabla\phi}\|_{",str,"}$"]);
            L1{3}      = join(["$\|\bar{\tau}_{f}^{\phantom{\nabla\phi}}\|_{",str,"}$"]);
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,NNZ,V1);
            %  > #2.
            M2         = repelem(":v",3);
            L2{1}      = join(["$\|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla}\phi}\|_{",str,"}$"]);
            L2{2}      = join(["$\|\bar{\tau}_{f^{\left(p\right)}}^{\nabla\phi}\|_{",str,"}$"]);
            L2{3}      = join(["$\|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla\phi}}\|_{",str,"}$"]);
            [L2,P2,Y2] = Fig_Tools_1D.Var_1(fig,M2,L2,NNZ,V2);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,[L1,L2],[P1,P2],NNZ,[Y1;Y2],[-1,1],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(edt,"Fig_V2_1_1",fig.Folder);
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_1_2(inp,obj,fig,j,n,edt)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            str   = Fig_Tools_1D.Set_str_1(j);
            p     = inp.ps.p+1;
            
            %  > Set variables to plot...
            for i = 1:size(obj,2)
                NNZ(i,1) = obj(i).m.nnz.At;
                V1 (i,1) = obj(i).e.a{n}.c.n_abs  (j);
                V2 (i,1) = obj(i).e.a{n}.t.n_abs.c(j);
                V3 (i,1) = obj(i).e.p{n}.c.n_abs  (j);
                V4 (i,1) = obj(i).e.p{n}.t.n_abs.c(j);
            end
            %  > Plot variables.
            %  > #1.
            M1         = repelem("--o",2);
            L1{1}      = join(["$\|e_{c}\|_{",str,"}$"]);
            L1{2}      = join(["$\|\bar{\tau}_{c}\|_{",str,"}$"]);
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,NNZ,[V1,V2]);
            %  > #2.
            M2         = repelem(":v",2);
            L2{1}      = join(["$\|e_{c^{\left(p\right)}}\|_{",str,"}$"]);
            L2{2}      = join(["$\|\bar{\tau}_{c^{\left(p\right)}}\|_{",str,"}$"]);
            [L2,P2,Y2] = Fig_Tools_1D.Var_1(fig,M2,L2,NNZ,[V3,V4]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,[L1,L2],[P1,P2],NNZ,[Y1;Y2],[-1,1],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(edt,"Fig_V2_1_2",fig.Folder);
        end
    end
end