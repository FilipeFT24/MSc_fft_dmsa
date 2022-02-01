classdef Fig_2_1D
    methods (Static)         
        %% > Wrap-up Fig_2 (1D).
        function [] = WrapUp_Fig_2_1D(Plot_2,Exp_2,Fig,msh,X,Norm)
            if Plot_2
                %  > Select...
                iF  = randperm(msh.f.NF,1);               
                %  > Figure.
                figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
                Fig_2_1D.Plot(msh,X,Norm);
                %  > Export as .pdf.
                if Exp_2
                    Fig_Tools.Export_PDF('Fig_2','../[Figures]/Fig_2');
                end
            end
        end
        
        %% > Auxiliary functions.
        function [] = Plot(msh,X,Norm)
            subplot(2,1,1);
            set(colorbar,'Visible','off'); 
            hold on;
            P1 = plot(msh.c.Xc,X.Error,'-ob','Linewidth',1.0,'MarkerFaceColor','b','MarkerSize',3.5);
            P2 = line([msh.c.Xc(1),msh.c.Xc(msh.c.NC)],[cell2mat(Norm.E(1)),cell2mat(Norm.E(1))],'Color','k','LineStyle','--','Linewidth',0.75);
            legend([P1,P2],["$\textrm{Absolute error}\left(\epsilon_{abs}\right)$","$|\!|\epsilon|\!|_{1}$",],'Interpreter','latex','Location','Northwest','FontSize',10);
            Fig_Tools_1D.ChangeLook_1(msh.f.Xv(1),msh.f.Xv(msh.f.NF),10);
        
            
            subplot(2,1,2);
            
            [ToPatch] = Fig_Tools_1D.ToPatch_Cell_Face(msh,0.05);
            hold on;
            for i = 1:msh.c.NC
                %  > NOTE: Add "'Linestyle','None'" to remove cell border.
                patch(ToPatch.Cell{i}(1,:),ToPatch.Cell{i}(2,:),X.Error(i));
            end
            c = Fig_Tools_1D.Colormap_cmocean('thermal');
            c.Label.String = '$\textrm{Absolute error}\left(\epsilon_{abs}\right)$';
            ax = gca; ax.YAxis.Visible = 'off'; 
            Fig_Tools_1D.ChangeLook_1(msh.f.Xv(1),msh.f.Xv(msh.f.NF),10);        
        end
    end
end      