classdef Fig_2_1D
    methods (Static)         
        %% > Wrap-up Fig_2 (1D).
        function [] = WrapUp_Fig_2_1D(Plot_2,Exp_2,Fig,msh,pde)
            if Plot_2               
                %  > Figure.
                figure(Fig); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot(msh,pde);
                %  > Export as .pdf.
                if Exp_2
                    Fig_Tools_1D.Export_PDF('Fig_2','../[Figures]/[1D]/Fig_2');
                end
            end
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [] = Plot(msh,pde)
            C  = linspecer(9,'qualitative');
            hold on;
            P1 = plot (msh.f.Xv,pde.en.f.f(:,3),'-^','Color',C(1,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            P2 = plot (msh.f.Xv,pde.en.f.f(:,1),'-v','Color',C(2,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            P3 = plot (msh.f.Xv,pde.en.f.f(:,2),'-o','Color',C(3,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            P4 = yline(pde.en.f.n(1)           ,'-' ,'Color',C(1,:),'Linewidth',1.0);
            P5 = yline(pde.en.f.n(3)           ,'-.','Color',C(1,:),'Linewidth',1.0);
            L  = Fig_2_1D.Set_Labels();
            set(colorbar,'visible','off');
            legend([P1,P2,P3,P4,P5],[L{1},L{2},L{3},L{4},L{5}],...
                'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',2);
            Fig_Tools_1D.ChangeLook_1D(true,true,msh.f.Xv,10,"$x$",L{6},20,12);
        end
        % >> 1.2. ---------------------------------------------------------
        function [L] = Set_Labels()
            L{1} = "$\epsilon_{f}^{\phi}$";
            L{2} = "$\bar{\epsilon}_{f}^{\phi}$";
            L{3} = "$\bar{\epsilon}_{f}^{\nabla\phi}$";
            L{4} = "$|\!|\epsilon_{f}^{\phi}|\!|_{1}$";
            L{5} = "$|\!|\epsilon_{f}^{\phi}|\!|_{\infty}$";
            L{6} = "$\textrm{Error magnitude}\left(\epsilon^{\phi}\right)$";
        end
    end
end      