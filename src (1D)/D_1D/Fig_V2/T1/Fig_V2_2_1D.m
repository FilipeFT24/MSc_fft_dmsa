classdef Fig_V2_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(V)
            %  > Auxiliary variables.
            exp  = 0;
            fig  = Fig_Tools_1D.Set_fig(1,exp);
            k    = 1;
            n    = 1;
            a    = 1:5;
            b    = 6:10;

            if ~exp
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_V2_2_1D.Plot_1_1(V(a),fig,k,n,0);
                subplot(1,2,2);
                Fig_V2_2_1D.Plot_1_1(V(b),fig,k,n,0);
            else
                plot = [0,1];
                %  > #1.
                if plot(1)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    Fig_V2_2_1D.Plot_1_1(V(a),fig,k,n,0);
                end
                %  > #2.
                if plot(2)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    Fig_V2_2_1D.Plot_1_1(V(b),fig,k,n,0);
                end
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1_1(V,fig,k,n,zoom)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            str   = Fig_Tools_1D.Set_str_1(k);
                        
            %  > Set variables to plot...
            for i = 1:size(V,2)
                for j = 1:size(V(i).obj,2)
                    h  (j,i) = V(i).msh(j).d.h;
                    nnz(j,i) = V(i).obj(j).m.nnz.At;
                    ec (j,i) = V(i).obj(j).e.a{n}.c.n_abs(k);
                    et (j,i) = V(i).obj(j).e.a{n}.t.n_abs.f(k,3);  
                end
                %  > Compute regression model...
                %  RM.ec{i} = C_1D.Set_RM(h(:,i),ec(:,i));
                %  RM.et{i} = C_1D.Set_RM(h(:,i),et(:,i));
            end
            %  > Plot variables.
            %  > #1.
            M1         = repelem("-o",size(nnz,2));
            L1{1}      = join(["$\|e_{c}^{\phantom{0}\left(2\right)}\|_{",str,"}$"]);
            L1{2}      = join(["$\|e_{c}^{\phantom{0}\left(4\right)}\|_{",str,"}$"]);
            L1{3}      = join(["$\|e_{c}^{\phantom{0}\left(6\right)}\|_{",str,"}$"]);
            L1{4}      = join(["$\|e_{c}^{\phantom{0}\left(8\right)}\|_{",str,"}$"]);
            L1{5}      = join(["$\|e_{c}^{\left(10\right)}\|_{",str,"}$"]);
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,nnz,ec);
            %  > #2.
            M2         = repelem(":v",size(nnz,2));
            L2{1}      = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{\phantom{0}\left(2\right)}\|_{",str,"}$"]);
            L2{2}      = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{\phantom{0}\left(4\right)}\|_{",str,"}$"]);
            L2{3}      = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{\phantom{0}\left(6\right)}\|_{",str,"}$"]);
            L2{4}      = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{\phantom{0}\left(8\right)}\|_{",str,"}$"]);
            L2{5}      = join(["$\|\bar{\tau}_{f^{\left(a\right)}}^{\left(10\right)}\|_{",str,"}$"]);
            [L2,P2,Y2] = Fig_Tools_1D.Var_1(fig,M2,L2,nnz,et);
            %  > Plot decay...
            %  p(1,:) = [200,2.0E-06];
            %  p(2,:) = [200,4.0E-08];
            %  p(3,:) = [200,5.0E-10];
            %  p(4,:) = [200,7.0E-12];
            %  p(5,:) = [200,1.0E-13];
            %  if fig.exp
            %      for i = 1:size(nnz,2)
            %          y (i) = 2.*i;
            %          X (i) = Fig_Tools_1D.Set_f(25,5.0E-5,y(i));
            %          P3{i} = fplot(@(x) X(i).*x.^-y(i),'Color',fig.C(i,:),'LineWidth',fig.LW./2,'Linestyle','--');
            %          text(p(i,1),p(i,2),join(['$$\mathcal{O}^{\left(',num2str(y(i)),'\right)}$$']),'FontSize',fig.FT_1,'Color',fig.C(i,:),'Interpreter','latex');
            %      end
            %  end
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,[L1,L2],[P1,P2],nnz,[Y1;Y2],[-1,4.5],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,"Fig_V2_2_1D");
        end
    end
end