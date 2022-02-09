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
            Xi = msh.f.Xv(1);
            Xf = msh.f.Xv(end);
            C  = linspecer(9,'qualitative');
            hold on;
            P1 = plot(msh.f.Xv,pde.e.f.f(:,1)           ,'-v','Color',C(1,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            P2 = plot(msh.f.Xv,pde.e.f.f(:,2)           ,'-o','Color',C(2,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            P3 = line([Xi,Xf],[pde.e.f.n(1,1),pde.e.f.n(1,1)],'Color',C(1,:),'Linewidth',1.0,'Linestyle','--');
            P4 = line([Xi,Xf],[pde.e.f.n(1,2),pde.e.f.n(1,2)],'Color',C(2,:),'Linewidth',1.0,'Linestyle','-.');
            L  = Fig_2_1D.Set_Labels();
            set(colorbar,'visible','off');
            legend([P1,P2,P3,P4],[L{1},L{2},L{3},L{4}],...
                'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',2);
            Fig_Tools_1D.ChangeLook_1D(true,true,msh.f.Xv,10,"$x$",L{5},20,12);
        end
        % >> 1.2. ---------------------------------------------------------
        function [L] = Set_Labels()
            L{1} = "$\bar{\epsilon}_{f}^{\phi}$";
            L{2} = "$\bar{\epsilon}_{f}^{\nabla\phi}$";
            L{3} = "$|\!|\bar{\epsilon}_{f}^{\phi}|\!|_{1}$";
            L{4} = "$|\!|\bar{\epsilon}_{f}^{\nabla\phi}|\!|_{1}$";
            L{5} = "$\textrm{Error magnitude}\left(\epsilon^{\phi}\right)$";
        end
    end
end      