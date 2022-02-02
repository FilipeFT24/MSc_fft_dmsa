classdef Fig_2_1D
    methods (Static)         
        %% > Wrap-up Fig_2 (1D).
        function [] = WrapUp_Fig_2_1D(Plot_2,Exp_2,Fig,msh,pde)
            if Plot_2
                %  > Select...
                iF  = randperm(msh.f.NF,1);               
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
            subplot(2,1,1);
            hold on;
            P1    = plot(msh.c.Xc,pde.E.EA,'-ob','Linewidth',1.0,'MarkerFaceColor','b','MarkerSize',2.0);
            P2    = line([msh.c.Xc(1),msh.c.Xc(msh.c.NC)],[pde.E.EN(1),pde.E.EN(1)],'Color','r','LineStyle','--','Linewidth',0.50);
            str_1 = "$\textrm{Absolute error}\left(\epsilon_{abs}\right)$";
            str_2 = '$|\!|\epsilon|\!|_{1}$';
            legend(P2,str_2,'Interpreter','latex','Location','Northwest','FontSize',10);
            %  > Colormap.
            set(colorbar,'Visible','off');
            %  > Axis.
            ylabel(str_1,'FontSize',20,'Interpreter','latex');
            Fig_Tools_1D.ChangeLook_1D(true,msh.f.Xv,10,12);            
            
            % >> Subplot 2.
            subplot(2,1,2);
            c_xy = Fig_Tools_1D.ToPatch_Cell_Face(msh,0.05);
            hold on;
            for i = 1:msh.c.NC
                %  > NOTE: Add "'Linestyle','None'" to remove cell border.
                patch(c_xy{i}(1,:),c_xy{i}(2,:),pde.E.EA(i),'Linestyle','None');
            end
            %  > Colormap.
            str_3 = '$\textrm{Absolute error}\left(\epsilon_{abs}\right)$';
            str_4 = 'thermal';
            c     = Fig_Tools_1D.Colormap_style(str_3,str_4,[0,roundn(max(pde.E.EA),ceil(log10(max(pde.E.EA))))],5,12);
            %  > Axis.
            ax = gca; ax.YAxis.Visible = 'off';
            Fig_Tools_1D.ChangeLook_1D(false,msh.f.Xv,10,12);        
        end
    end
end      