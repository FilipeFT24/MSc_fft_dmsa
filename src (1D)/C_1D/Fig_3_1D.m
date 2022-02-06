classdef Fig_3_1D
    methods (Static)         
        %% > Wrap-up Fig_3 (1D).
        function [] = WrapUp_Fig_3_1D(Plot_3,Exp_3,Fig,msh,pde)
            if Plot_3               
                %  > Figure.
                figure(Fig); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_3_1D.Plot(msh,pde);
                %  > Export as .pdf.
                if Exp_3
                    Fig_Tools_1D.Export_PDF('Fig_3','../[Figures]/[1D]/Fig_3');
                end
            end
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [] = Plot(msh,pde)
            % >> 1.1.1. ---------------------------------------------------
            subplot(2,1,1);
            C  = linspecer(9,'qualitative');
            hold on;
            P1 = plot (msh.c.Xc,pde.en.c.c(:,1),'-s','Color',C(1,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            P2 = yline(pde.en.c.n(1)           ,'-' ,'Color',C(1,:),'Linewidth',1.0);
            L  = Fig_3_1D.Set_Labels();
            set(colorbar,'visible','off');
            legend([P1,P2],[L{1},L{2}],...
                'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',2);
            Fig_Tools_1D.ChangeLook_1D(true,true,msh.f.Xv,10,"$x$",L{3},20,12);
            
            % >> 1.1.2. ---------------------------------------------------
            subplot(2,1,2);
            c_xy = Fig_Tools_1D.ToPatch(msh,0.01);
            hold on;
            for i = 1:msh.c.NC
                patch(c_xy{i}(1,:),c_xy{i}(2,:),pde.en.c.c(i),'Linestyle','None');
            end
            %  > Colormap.
            c  = Fig_Tools_1D.Colormap_style(L{3},'thermal',12);
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(false,false,msh.f.Xv,10,"$x$","$y$",20,12);    
            ax = gca; ax.YAxis.Visible = 'off';
        end
        % >> 1.2. ---------------------------------------------------------
        function [L] = Set_Labels()
            L{1} = "$\epsilon_{c}^{\phi}$";
            L{2} = "$|\!|\epsilon_{c}^{\phi}|\!|_{1}$";
            L{3} = "$\textrm{Error magnitude}\left(\epsilon^{\phi}\right)$";
        end
    end
end      