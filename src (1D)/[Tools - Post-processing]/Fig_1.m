classdef Fig_1
    methods (Static)         
        %% > Wrap up Fig_1.
        function [] = WrapUp_Fig_1(msh)
            if strcmpi(obj.msh_Type,'Non-uniform')
                figure(1); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
                %  > Subplot 1.
                Fig_1.Subplot_1(msh);
                %  > Subplot 2.
                Fig_1.Subplot_2(msh);
            end
        end
        
        %% > Auxiliary functions.
        %  > Subplot 1.
        function [] = Subplot_1(msh)
            %  > Uniform mesh cell spacement.
            dX = (obj.Xv_f-obj.Xv_i)./msh.NC;
            %  > Cell volume ratio(r_1).
            for i = 1:msh.NC
                r_1(i) = msh.Vol(i)./dX;
            end
            %  > Cell(i+1)/Cell(i) volume ratio(r_2).
            for i = 1:msh.NC-1
                x_2(i) = msh.Xv(i+1);
                r_2(i) = msh.Vol(i+1)./msh.Vol(i);
            end
            subplot(2,1,1);
            hold on;
            P1 = plot(msh.Xc,r_1,'-or','MarkerFaceColor','r','MarkerSize',3.5);
            P2 = plot(x_2,r_2   ,'-ob','MarkerFaceColor','b','MarkerSize',3.5);
            legend([P1,P2],["$\Delta\mathrm{x}/\Delta\mathrm{x}_{0}$","$\Delta\mathrm{x_{_{RHS}}}/\Delta\mathrm{x_{_{LHS}}}$",],'Interpreter','latex','Location','Northwest','FontSize',10);
            Fig_Tools.ChangeLook_1(msh.Xv(1),msh.Xv(msh.NV),10);
        end
        %  > Subplot 2.
        function [] = Subplot_2(msh)
            %  > To patch...
            [ToPatch] = Fig_Tools.ToPatch_Cell_Face(msh,0.05);
            subplot(2,1,2);
            hold on;
            for i = 1:msh.NC
                %  > NOTE: Add "'Linestyle','None'" to remove cell border.
                patch(ToPatch.Cell{i}(1,:),ToPatch.Cell{i}(2,:),'w');
            end
            ax = gca; ax.YAxis.Visible = 'off'; 
            Fig_Tools.ChangeLook_1(msh.Xv(1),msh.Xv(msh.NV),10); 
        end
    end
end         