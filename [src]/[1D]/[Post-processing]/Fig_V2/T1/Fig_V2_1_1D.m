classdef Fig_V2_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(V)
            %  > Auxiliary variables.
            exp  = 1;
            fig  = Fig_Tools_1D.Set_fig(1,exp);
            k    = 1;
            n    = 1;
            plot = [1,0];

            if ~exp
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_V2_1_1D.Plot_1_1(V(1:2),fig,k,n,0);
                subplot(1,2,2);
                Fig_V2_1_1D.Plot_1_1(V(3:4),fig,k,n,0);
            else
                %  > #1.
                if plot(1)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    Fig_V2_1_1D.Plot_1_1(V(1:2),fig,k,n,0);
                end
                %  > #2.
                if plot(2)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    Fig_V2_1_1D.Plot_1_1(V(3:4),fig,k,n,0);
                end
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1_1(V,fig,k,n,zoom)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            str_1 = ["1/n-1","1/n"];
            str_2 = Fig_Tools_1D.Set_str_1(k);
                        
            %  > Set variables to plot...
            for i = 1:size(V,2)
                for j = 1:size(V(i).obj,2)
                    h     (j,i) = V(i).msh(j).d.h;
                    nnz   (j,i) = V(i).obj(j).m.nnz.At;
                    ec    (j,i) = V(i).obj(j).e.a{n}.c.n_abs(k);
                    et {i}(j,:) = V(i).obj(j).e.a{n}.t.n_abs.f(k,:);  
                end
            end
            NNZ = nnz(:,1);
            %  > Plot variables.
            %  > #1.
            l          = 1;
            M1         = ["-o","-o","-o","-d"];
            L1{1}      = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{\phantom{\nabla}\phi_{\left(",str_1(l),"\right)}}\|_{",str_2,"}$"]);
            L1{2}      = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{\nabla\phi_{\left(",str_1(l),"\right)}}\|_{",str_2,"}$"]);
            L1{3}      = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{\phantom{\nabla\phi}_{\left(",str_1(l),"\right)}}\|_{",str_2,"}$"]);
            L1{4}      = join(["$\|e_{c^{\left(a\right)}}^{\phantom{\nabla\phi}_{\left(",str_1(l),"\right)}}\|_{",str_2,"}$"]);
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,NNZ,[et{l},ec(:,l)]);
            %  > #2.
            l          = 2;
            M2         = [":v",":v",":v",":d"];
            L2{1}      = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{\phantom{\nabla}\phi_{\left(",str_1(l),"\right)}}\|_{",str_2,"}$"]);
            L2{2}      = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{\nabla\phi_{\left(",str_1(l),"\right)}}\|_{",str_2,"}$"]);
            L2{3}      = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{\phantom{\nabla\phi}_{\left(",str_1(l),"\right)}}\|_{",str_2,"}$"]);
            L2{4}      = join(["$\|e_{c^{\left(a\right)}}^{\phantom{\nabla\phi}_{\left(",str_1(l),"\right)}}\|_{",str_2,"}$"]);
            [L2,P2,Y2] = Fig_Tools_1D.Var_1(fig,M2,L2,NNZ,[et{2},ec(:,2)]);
            %  > Plot decay...
            %  if fig.exp
            %      %  > #1.
            %      c      = min(max(ec))./15;
            %      p      = [375,4.5E-6];
            %      %  > #2.
            %      c      = min(max(ec))./80;
            %      p      = [375,7.5E-9];
            %      X3     = Fig_Tools_1D.Set_f(NNZ(1),c,2);
            %      P3     = fplot(@(x) X3.*x.^-2,'--k','LineWidth',fig.LW./2);
            %      text(p(1),p(2),join(['$$\mathcal{O}^{\left(',num2str(2),'\right)}$$']),'FontSize',fig.FT_1,'Interpreter','latex');
            %  end
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,[L1,L2],[P1,P2],[nnz,nnz],[Y1;Y2],[-1,0],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,"Fig_V2_1_1D");
        end
    end
end