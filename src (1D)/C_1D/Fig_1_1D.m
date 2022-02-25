classdef Fig_1_1D
    methods (Static)
        %% > Wrap-up Fig_1 (1D).
        function [] = WrapUp_Fig_1_1D(msh,pde)
            %  > Plot/Export.
            P_1 = true;
            P_2 = true;
            P_3 = true;
            Exp = false;
            F_1 = "FN_1";
            F_2 = "FN_2";
            F_2 = "FN_3";
            D_1 = "../[Figures]/[1D]/Fig_2";
            D_2 = "../[Figures]/[1D]/Fig_2";
            D_3 = "../[Figures]/[1D]/Fig_2";
            %  > Properties.
            Fig = [1,2,3];
            fig = Fig_1_1D.Set_fig(Exp);
            l_1 = true;
            l_2 = true;
            
            if ~Exp
                figure(Fig(1)); set(gcf,'Units','pixels','Position',[150,100,1250,600]);
                subplot(1,2,1);
                Fig_1_1D.Plot_1(fig,l_1,msh,pde);
                subplot(1,2,2);
                Fig_1_1D.Plot_2(fig,l_2,msh,pde);
            else
                if P_1
                    figure(Fig(1)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                    Fig_1_1D.Plot_1(fig,l_1,msh,pde);
                    Fig_Tools_1D.Export_PDF(F_1,D_1);
                end
                if P_2
                    figure(Fig(2)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                    Fig_1_1D.Plot_2(fig,l_2,msh,pde);
                    Fig_Tools_1D.Export_PDF(F_2,D_2);
                end
                if P_3
                    figure(Fig(3)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                    Fig_1_1D.Plot_3(fig,l_3,msh,pde);
                    Fig_Tools_1D.Export_PDF(F_3,D_3);
                end
            end
        end
        
        %% > 1. -----------------------------------------------------------
        function [] = Plot_1(fig,l_1,msh,pde)
            %  > Auxiliary variables.
            C  = linspecer(9,'qualitative');
            L  = Fig_1_1D.Set_Labels_1();
            m  = size(pde.e.t.f,2);
            Xv = [msh.f.Xv(1),msh.f.Xv(msh.f.NF)];
            
            %  > Plot.
            hold on;
            for j = 1:m+1
                if j ~= m+1
                    P{j}     = plot(msh.f.Xv,pde.e.t.f(:,j)          ,':o','Color',C(j,:),'LineWidth',fig.LW_1,'MarkerFaceColor',C(j,:),'MarkerSize',fig.MS_1);
                    P{j+m+1} = line(Xv,[pde.e.t.n.f(1,j),pde.e.t.n.f(1,j)],'Color',C(j,:),'Linewidth',fig.LW_2,'Linestyle','-.');
                else
                    P{j}     = plot(msh.c.Xc,pde.e.t.c(:,1)          ,':o','Color',C(j,:),'LineWidth',fig.LW_1,'MarkerFaceColor',C(j,:),'MarkerSize',fig.MS_1);
                    P{j+m+1} = line(Xv,[pde.e.t.n.c(1,1),pde.e.t.n.c(1,1)],'Color',C(j,:),'Linewidth',fig.LW_2,'Linestyle','-.');
                end
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',fig.FT_1,'NumColumns',2);
            %  > Axis.
            if l_1
                set(gca,'YScale','log');
            end
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{7},fig.FT_2,fig.FT_3);
        end
        
        %% > 2. -----------------------------------------------------------
        function [] = Plot_2(fig,l_2,msh,pde)
            %  > Auxiliary variables.
            C  = linspecer(9,'qualitative');
            L  = Fig_1_1D.Set_Labels_2();
            m  = size(pde.e.f.f,2);
            Xv = [msh.f.Xv(1),msh.f.Xv(msh.f.NF)];
            
            %  > Plot.
            hold on;
            for j = 1:m+1
                if j ~= m+1
                    P{j}     = plot(msh.f.Xv,pde.e.f.f(:,j)      ,':o','Color',C(j,:),'LineWidth',fig.LW_1,'MarkerFaceColor',C(j,:),'MarkerSize',fig.MS_1);
                    P{j+m+1} = line(Xv,[pde.e.f.n(1,j),pde.e.f.n(1,j)],'Color',C(j,:),'Linewidth',fig.LW_2,'Linestyle','-.');
                else
                    P{j}     = plot(msh.c.Xc,pde.e.c.c(:,1)      ,':o','Color',C(j,:),'LineWidth',fig.LW_1,'MarkerFaceColor',C(j,:),'MarkerSize',fig.MS_1);
                    P{j+m+1} = line(Xv,[pde.e.c.n(1,1),pde.e.c.n(1,1)],'Color',C(j,:),'Linewidth',fig.LW_2,'Linestyle','-.');
                end
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',fig.FT_1,'NumColumns',2);
            %  > Axis.
            if l_2
                set(gca,'YScale','log');
            end
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{7},fig.FT_2,fig.FT_3);
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        function [L] = Set_Labels_1()
            L{1} = "$|\tau_{f}^{\phantom{\nabla}\phi}|$";
            L{2} = "$|\tau_{f}^{\nabla\phi}|$";
            L{3} = "$|\tau_{c}^{\phantom{\nabla\phi}}|$";
            L{4} = "$|\!|\tau_{f}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{5} = "$|\!|\tau_{f}^{\nabla\phi}|\!|_{1}$";
            L{6} = "$|\!|\tau_{c}^{\phantom{\nabla\phi}}|\!|_{1}$";
            L{7} = "$\textrm{Error magnitude}, |\tau|$";
        end
        % >> 3.2. ---------------------------------------------------------
        function [L] = Set_Labels_2()
            L{1} = "$|e_{f}^{\phantom{\nabla}\phi}|$";
            L{2} = "$|e_{f}^{\nabla\phi}|$";
            L{3} = "$|e_{c}^{\phantom{\nabla\phi}}|$";
            L{4} = "$|\!|e_{f}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{5} = "$|\!|e_{f}^{\nabla\phi}|\!|_{1}$";
            L{6} = "$|\!|e_{c}^{\phantom{\nabla\phi}}|\!|_{1}$";
            L{7} = "$\textrm{Error magnitude}, |e|$";
        end
        % >> 3.3. ---------------------------------------------------------
        function [fig] = Set_fig(Exp)
            if ~Exp
                fig.LW_1 = 2.0;
                fig.LW_2 = 1.5;
                fig.MS_1 = 3.0;
                fig.FT_1 = 12.5;
                fig.FT_2 = 12.5;
                fig.FT_3 = 12.5;
            else
                fig.LW_1 = 3.0;
                fig.LW_2 = 2.5;
                fig.MS_1 = 5.0;
                fig.FT_1 = 25.0;
                fig.FT_2 = 30.0;
                fig.FT_3 = 30.0;
            end
        end
    end
end