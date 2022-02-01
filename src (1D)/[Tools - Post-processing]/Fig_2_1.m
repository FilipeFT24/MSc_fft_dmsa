classdef Fig_2_1
    methods (Static)         
        %% > Wrap up Fig_2_1.
        function [] = WrapUp_Fig_2_1(msh,X,Norm)
            figure(2); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
            if strcmpi(obj.Sim_Type,'Explicit')
                %  > Subplot 1.
                Fig_2_1.Subplot_1_EX(msh,X);
                %  > Subplot 2.
                Fig_2_1.Subplot_2_EX(msh,X);
                %  > Subplot 3.
                Fig_2_1.Subplot_3_EX(msh,X,Norm);
                %  > Subplot 4.
                Fig_2_1.Subplot_4_EX(msh,X);
            elseif strcmpi(obj.Sim_Type,'Implicit') || strcmpi(obj.Sim_Type,'DC')
                %  > Subplot 1.
                Fig_2_1.Subplot_1(msh,X);
                %  > Subplot 2.
                Fig_2_1.Subplot_2(msh,X);
                %  > Subplot 3.
                Fig_2_1.Subplot_3(msh,X,Norm);
                %  > Subplot 4.
                Fig_2_1.Subplot_4(msh,X);
            end
        end
        
        %% > Auxiliary functions.
        % >> Explicit flux reconstruction. 
        %  > Subplot 1.
        function [] = Subplot_1_EX(msh,X)
            subplot(4,1,1);
            set(colorbar,'Visible','off','Location','EastOutside'); 
            hold on;
            plot(msh.Xc,X.Phi,'-k','Linewidth',1.0);
            plot(msh.Xc,X.Phi,'or','MarkerFaceColor','r','MarkerSize',3.5);
            ylabel("$\phi$",'Interpreter','latex','FontSize',10);
            Fig_Tools.ChangeLook_1(msh.Xv(1),msh.Xv(msh.NV),10);
        end
        %  > Subplot 2.
        function [] = Subplot_2_EX(msh,X)
            %  > To patch...
            [ToPatch] = Fig_Tools.ToPatch_Cell_Face(msh,0.05);
            subplot(4,1,2);
            hold on;
            for i = 1:msh.NC
                %  > NOTE: Add "'Linestyle','None'" to remove cell border.
                patch(ToPatch.Cell{i}(1,:),ToPatch.Cell{i}(2,:),X.Phi(i));
            end
            c = Fig_Tools.Colormap_cmocean();
            c.Label.String = '$\phi$';
            ax = gca; ax.YAxis.Visible = 'off'; 
            Fig_Tools.ChangeLook_1(msh.Xv(1),msh.Xv(msh.NV),10);             
        end
        %  > Subplot 3.
        function [] = Subplot_3_EX(msh,X,Norm)
            subplot(4,1,3);
            set(colorbar,'Visible','off'); 
            hold on;
            P1 = plot(msh.Xc,X.Error,'-ob','Linewidth',1.0,'MarkerFaceColor','b','MarkerSize',3.5);
            P2 = line([msh.Xc(1),msh.Xc(msh.NC)],[cell2mat(Norm.E(1)),cell2mat(Norm.E(1))],'Color','k','LineStyle','--','Linewidth',0.75);
            legend([P1,P2],["$\textrm{Absolute error}\left(\epsilon_{abs}\right)$","$|\!|\epsilon|\!|_{1}$",],'Interpreter','latex','Location','Northwest','FontSize',10);
            Fig_Tools.ChangeLook_1(msh.Xv(1),msh.Xv(msh.NV),10);
        end
        %  > Subplot 4.
        function [] = Subplot_4_EX(msh,X)
            %  > To patch...
            [ToPatch] = Fig_Tools.ToPatch_Cell_Face(msh,0.05);
            subplot(4,1,4);
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
        
        % >> Implicit/DC approach flux reconstruction. 
        %  > Subplot 1.
        function [] = Subplot_1(msh,X)
            subplot(4,1,1);
            set(colorbar,'Visible','off'); 
            hold on;
            P1 = plot(msh.Xc,X.Phi    ,'-k','Linewidth',1.0);
            P2 = plot(msh.Xc,X.Phi_PDE,'or','MarkerFaceColor','r','MarkerSize',3.5);
            legend([P1,P2],["$\phi$","$\phi_{_{\mathrm{PDE}}}$",],'Interpreter','latex','Location','Northwest','FontSize',10);
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
                patch(ToPatch.Cell{i}(1,:),ToPatch.Cell{i}(2,:),X.Phi_PDE(i));
            end
            c = Fig_Tools.Colormap_cmocean();
            c.Label.String = '$\phi_{_{\mathrm{PDE}}}$'; 
            ax = gca; ax.YAxis.Visible = 'off';
            Fig_Tools.ChangeLook_1(msh.Xv(1),msh.Xv(msh.NV),10);            
        end
        %  > Subplot 3.
        function [] = Subplot_3(msh,X,Norm)
            subplot(4,1,3);
            set(colorbar,'Visible','off'); 
            hold on;
            P1 = plot(msh.Xc,X.Error,'-ob','Linewidth',1.0,'MarkerFaceColor','b','MarkerSize',3.5);
            P2 = line([msh.Xc(1),msh.Xc(msh.NC)],[cell2mat(Norm.E(1)),cell2mat(Norm.E(1))],'Color','k','LineStyle','--','Linewidth',0.75);
            legend([P1,P2],["$\textrm{Absolute error}\left(\epsilon_{abs}\right)$","$|\!|\epsilon|\!|_{1}$",],'Interpreter','latex','Location','Northwest','FontSize',10);
            Fig_Tools.ChangeLook_1(msh.Xv(1),msh.Xv(msh.NV),10);
        end
        %  > Subplot 4.
        function [] = Subplot_4(msh,X)
            %  > To patch...
            [ToPatch] = Fig_Tools.ToPatch_Cell_Face(msh,0.05);
            subplot(4,1,4);
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
    end
end      