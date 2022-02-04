classdef Fig_2_1D
    methods (Static)         
        %% > Wrap-up Fig_2 (1D).
        function [] = WrapUp_Fig_2_1D(Plot_2,Exp_2,Fig,msh,pde)
            if Plot_2               
                %  > Figure.
                figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
                Fig_2_1D.Plot(msh,pde);
                %  > Export as .pdf.
                if Exp_2
                    Fig_Tools_1D.Export_PDF('Fig_2','../[Figures]/[1D]/Fig_2');
                end
            end
        end
        
        %% > Auxiliary functions.
        function [] = Plot(msh,pde)
            % >> Subplot 1.
            %  > Auxiliary array.
            for i = 1:msh.f.NF
                en_f(i,1) = pde.en.f{i}(1);
                en_f(i,2) = pde.en.f{i}(2);
            end
            
            subplot(2,1,1);
            hold on;
            C     = linspecer(3,'qualitative');
            %  > Cell(s).
            P1    = plot(msh.c.Xc,pde.en.c.c,'-s','Color',C(1,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            %  > Face(s).
            P2    = plot(msh.f.Xv,en_f(:,1) ,'-o','Color',C(2,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            P3    = plot(msh.f.Xv,en_f(:,2) ,'-^','Color',C(3,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            %  > Mean (cell) error.
            P4    = yline(pde.en.c.n(1)     ,'Color','k','LineStyle','--','Linewidth',0.5);
            %  > Label(s).
            str_1 = "$\epsilon_{c}^{\phi}$";
            str_2 = "$\epsilon_{f}^{\phi}$";
            str_3 = "$\epsilon_{f}^{\nabla\phi}$";
            str_4 = "$|\!|\epsilon_{c}^{\phi}|\!|_{1}$";
            str_5 = "$\textrm{Absolute error}\left(\epsilon_{abs}\right)$";
            legend([P1,P2,P3,P4],[str_1,str_2,str_3,str_4],'Interpreter','latex','Location','Northeastoutside','FontSize',10);
            %  > Axis.
            ylabel(str_5,'FontSize',20,'Interpreter','latex');
            Fig_Tools_1D.ChangeLook_1D(true,msh.f.Xv,10,12);            
            
            % >> Subplot 2.
            subplot(2,1,2);
            c_xy = Fig_Tools_1D.ToPatch(msh,0.05);
            hold on;
            for i = 1:msh.c.NC
                %  > NOTE: Add "'Linestyle','None'" to remove cell border.
                patch(c_xy{i}(1,:),c_xy{i}(2,:),pde.en.c.c(i),'Linestyle','None');
            end
            %  > Colormap.
            str_2 = '$\textrm{Absolute error}\left(\epsilon_{abs}\right)$';
            str_4 = 'thermal';
            c     = Fig_Tools_1D.Colormap_style(str_2,str_4,[0,roundn(max(pde.en.c.c),ceil(log10(max(pde.en.c.c))))],5,12);
            %  > Axis.
            ax = gca; ax.YAxis.Visible = 'off';
            Fig_Tools_1D.ChangeLook_1D(false,msh.f.Xv,10,12);        
        end
    end
end      