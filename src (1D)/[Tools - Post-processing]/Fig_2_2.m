classdef Fig_2_2
    methods (Static)         
        %% > Wrap up Fig_2_2.
        function [] = WrapUp_Fig_2_2(msh,X,Norm)
            figure(1); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
            %  > Subplot 1.
            Fig_2_2.Subplot_1(msh{1,1},X{1,1},Norm{1,1});
            %  > Subplot 2.
            Fig_2_2.Subplot_2(msh{1,1},X{1,1});
            %  > Subplot 3.
            Fig_2_2.Subplot_3(msh,X);
            %  > Subplot 4.
            Fig_2_2.Subplot_4(msh,X);
        end
        
        %% > Auxiliary functions.
        %  > Subplot 1.
        function [] = Subplot_1(msh,X,Norm)
            subplot(4,1,1);
            set(colorbar,'Visible','off'); 
            hold on;
            plot(msh.Xc,X.Error,'-ob','Linewidth',1.0,'MarkerFaceColor','b','MarkerSize',3.5);
            P1 = line([msh.Xc(1),msh.Xc(msh.NC)],[cell2mat(Norm.E(1)),cell2mat(Norm.E(1))],'Color','k','LineStyle','--','Linewidth',0.75);
            P2 = line([msh.Xc(1),msh.Xc(msh.NC)],[obj.K_ref.*cell2mat(Norm.E(1)),obj.K_ref.*cell2mat(Norm.E(1))],'Color','r','LineStyle','--','Linewidth',0.75);
            ylabel("$\textrm{Absolute error}\left(\epsilon_{abs}\right)$",'Interpreter','latex','FontSize',10);
            legend([P1,P2],["$|\!|\epsilon|\!|_{1}$","$|\!|\epsilon|\!|_{1}^*$"],'Interpreter','latex','Location','Northwest','FontSize',10);
            Fig_Tools.ChangeLook_1(msh.Xv(1),msh.Xv(msh.NV),10);
        end
        %  > Subplot 2.
        function [] = Subplot_2(msh,X)
            %  > To patch...
            [ToPatch] = Fig_Tools.ToPatch_Cell_Face(msh,0.05);
            
            subplot(4,1,2);
            hold on;
            for i = 1:msh.NC
                %  > NOTE: Add "'Linestyle','None'" to remove cell border.
                patch(ToPatch.Cell{i}(1,:),ToPatch.Cell{i}(2,:),X.Error(i));
            end
            c = Fig_Tools.Colormap_cmocean();
            c.Label.String = '$\textrm{Absolute error}\left(\epsilon_{abs}\right)$';
            ax = gca; ax.YAxis.Visible = 'off'; 
            Fig_Tools.ChangeLook_1(msh.Xv(1),msh.Xv(msh.NV),10); 
        end
        %  > Subplot 3.
        function [] = Subplot_3(msh,X)
            subplot(4,1,3);
            C = linspecer(size(msh,2),'qualitative');
            set(colorbar,'Visible','off'); 
            hold on;
            for i = 1:size(msh,2)
                P{i} = plot(msh{1,i}.Xc,X{1,i}.Error,'-^','color',C(i,:),'Linewidth',1.5,'MarkerFaceColor',C(i,:),'MarkerSize',1.5);
                legendName(i) = convertCharsToStrings(num2str(i));
            end
            ylabel("$\textrm{Absolute error}\left(\epsilon_{abs}\right)$",'Interpreter','latex','FontSize',10);
            legend(legendName,'Interpreter','latex','Location','Northeastoutside','FontSize',10);
            l = size(msh,2);
            Fig_Tools.ChangeLook_1(msh{1,l}.Xv(1),msh{1,l}.Xv(msh{1,l}.NV),10);
        end
        %  > Subplot 4.
        function [] = Subplot_4(msh,X)
            %  > Select.
            l = size(msh,2);
            %  > To patch...
            [ToPatch] = Fig_Tools.ToPatch_Cell_Face(msh{1,l},0.05);
                       
            subplot(4,1,4);
            hold on;
            for i = 1:msh{1,l}.NC
                %  > NOTE: Add "'Linestyle','None'" to remove cell border.
                patch(ToPatch.Cell{i}(1,:),ToPatch.Cell{i}(2,:),X{1,l}.Error(i));
            end
            c = Fig_Tools.Colormap_cmocean();
            c.Label.String = '$\textrm{Absolute error}\left(\epsilon_{abs}\right)$';
            ax = gca; ax.YAxis.Visible = 'off'; 
            Fig_Tools.ChangeLook_1(msh{1,l}.Xv(1),msh{1,l}.Xv(msh{1,l}.NV),10); 
        end
    end
end  