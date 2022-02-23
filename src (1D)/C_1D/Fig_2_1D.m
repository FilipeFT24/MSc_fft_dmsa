classdef Fig_2_1D
    methods (Static)
        %% > Wrap-up Fig_2 (1D).
        function [] = WrapUp_Fig_2_1D(msh,pde)
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
            fig = Fig_2_1D.Set_fig(Exp);
            l_1 = false;
            l_2 = false;
            
            if ~Exp
                figure(Fig(1)); set(gcf,'Units','pixels','Position',[150,100,1250,600]);
                subplot(1,2,1);
                Fig_2_1D.Plot_1(fig,l_1,msh,pde);
                subplot(1,2,2);
                Fig_2_1D.Plot_2(fig,l_2,msh,pde);
            else
                if P_1
                    figure(Fig(1)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                    Fig_2_1D.Plot_1(fig,l_1,msh,pde);
                    Fig_Tools_1D.Export_PDF(F_1,D_1);
                end
                if P_2
                    figure(Fig(2)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                    Fig_2_1D.Plot_2(fig,l_2,msh,pde);
                    Fig_Tools_1D.Export_PDF(F_2,D_2);
                end
                if P_3
                    figure(Fig(3)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                    Fig_2_1D.Plot_3(fig,l_3,msh,pde);
                    Fig_Tools_1D.Export_PDF(F_3,D_3);
                end
            end
        end
        
        %% > 1. -----------------------------------------------------------
        function [] = Plot_1(f_1,l_1,msh,pde)
            %  > Auxiliary variables.
            C  = linspecer(9,'qualitative');
            L  = Fig_2_1D.Set_Labels_1();
            m  = size(pde.e.t.f,2);
            Xv = [msh.f.Xv(1),msh.f.Xv(msh.f.NF)];
            
            %  > Plot.
            hold on;
            for j = 1:m
                P{j}   = plot(msh.f.Xv,pde.e.t.f(:,j)          ,':o','Color',C(j,:),'LineWidth',f_1.LW_1,'MarkerFaceColor',C(j,:),'MarkerSize',f_1.MS_1);
                P{j+m} = line(Xv,[pde.e.t.n.f(1,j),pde.e.t.n.f(1,j)],'Color',C(j,:),'Linewidth',f_1.LW_2,'Linestyle','-.');
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',f_1.FT_1,'NumColumns',2);
            %  > Axis.
            if l_1
                set(gca,'YScale','log');
            end
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{8},f_1.FT_2,f_1.FT_3);
        end
        
        %% > 2. -----------------------------------------------------------
        function [] = Plot_2(f_2,l_2,msh,pde)
            %  > Auxiliary variables.
            C  = linspecer(9,'qualitative');
            L  = Fig_2_1D.Set_Labels_2();
            Xv = [msh.f.Xv(1),msh.f.Xv(msh.f.NF)];
            
            %  > Plot.
            hold on;
            P{1} = plot(msh.c.Xc,pde.e.c.c(:,1)      ,':o','Color',C(1,:),'LineWidth',f_2.LW_1,'MarkerFaceColor',C(1,:),'MarkerSize',f_2.MS_1);
            P{2} = line(Xv,[pde.e.c.n(1,1),pde.e.c.n(1,1)],'Color',C(1,:),'Linewidth',f_2.LW_2,'Linestyle','-.');
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',f_2.FT_1,'NumColumns',2);
            %  > Axis.
            if l_2
                set(gca,'YScale','log');
            end
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{3},f_2.FT_2,f_2.FT_3);
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        function [L] = Set_Labels_1()
            L{1} = "$|\tau_{f}^{\phantom{\nabla}\phi}|$";
            L{2} = "$|\tau_{f}^{\nabla\phi}|$";
            L{3} = "$|\!|\tau_{f}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{4} = "$|\!|\tau_{f}^{\nabla\phi}|\!|_{1}$";
            L{5} = "$|\!|\bar\tau_{f}^{\phantom{\nabla\phi}}|\!|_{1}$";
            L{6} = "$|\tau_{c}|$";
            L{7} = "$|\!|\tau_{c}|\!|_{1}$";
            L{8} = "$\textrm{Error magnitude}, |\tau^{\phi}|$";
        end
        % >> 3.2. ---------------------------------------------------------
        function [L] = Set_Labels_2()
            L{1} = "$|e_{c}^{\phantom{\nabla}\phi}|$";
            L{2} = "$|\!|e_{c}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{3} = "$\textrm{Error magnitude}, |e^{\phi}|$";
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