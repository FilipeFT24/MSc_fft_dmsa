classdef Fig_4
    methods (Static)         
        %% > Wrap up Fig_4.
        function [] = WrapUp_Fig_4(msh,X1,X2)
            figure(2); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
            %  > Subplot 1.
            Fig_4.Subplot_1(msh,X1,X2);
            %  > Subplot 2/3/4.
            Fig_4.Subplot_2(msh,X1,X2.Richardson.Phi_L,2);
            Fig_4.Subplot_2(msh,X1,X2.Richardson.Phi_U,3);
            Fig_4.Subplot_2(msh,X1,X2.Richardson.Phi_C,4);
        end
        
        %% > Auxiliary functions. 
        % >> Subplot 1.
        function [] = Subplot_1(msh,X1,X2)
            subplot(4,1,1);
            C = linspecer(3,'qualitative');
            hold on;
            P1 = plot(msh{1}.Xc,abs(X1{1}.Phi-X2.Richardson.Phi_L),'-^','color',C(1,:),'Linewidth',1.5,'MarkerFaceColor',C(1,:),'MarkerSize',3.5);
            P2 = plot(msh{1}.Xc,abs(X1{1}.Phi-X2.Richardson.Phi_U),'-^','color',C(2,:),'Linewidth',1.5,'MarkerFaceColor',C(2,:),'MarkerSize',3.5);
            P3 = plot(msh{1}.Xc,abs(X1{1}.Phi-X2.Richardson.Phi_C),'-^','color',C(3,:),'Linewidth',1.5,'MarkerFaceColor',C(3,:),'MarkerSize',3.5);
            set(colorbar,'Visible','off');
            legend([P1,P2,P3],...
                ["$|\!|e_{_{L}}|\!|$","$|\!|e_{_{U}}|\!|$","$|\!|e_{_{C}}|\!|$"],'Interpreter','latex','Location','Northwest','FontSize',10);
            Fig_Tools.ChangeLook(msh{1}.Xv(1),msh{1}.Xv(msh{1}.NV),10);
        end
        % >> Subplot 2/3/4.
        function [] = Subplot_2(msh,X1,M,i)
            %  > To patch...
            [ToPatch] = Fig_Tools.ToPatch_Cell_Face(msh{1},0.05);
            subplot(4,1,i);
            hold on;
            for i = 1:msh{1}.NC
                %  > NOTE: Add "'Linestyle','None'" to remove cell border.
                patch(ToPatch.Cell{i}(1,:),ToPatch.Cell{i}(2,:),abs(X1{1}.Phi(i)-M(i)));
            end
            Fig_Tools.Colormap_cmocean();
            Fig_Tools.ChangeLook(msh{1}.Xv(1),msh{1}.Xv(msh{1}.NV),10);
        end
    end
end
        
        
        
        
        
        
        
        
        