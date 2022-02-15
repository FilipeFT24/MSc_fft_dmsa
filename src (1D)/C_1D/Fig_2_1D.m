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
            %  > Auxiliary variables.
            Xvi = msh.f.Xv(1);
            Xvf = msh.f.Xv(end);
            C   = linspecer(9,'qualitative');
            L   = Fig_2_1D.Set_Labels();
            %  > 1.1.1. ---------------------------------------------------
            subplot(2,1,1);
            Fig_2_1D.SubPlot_1(msh,pde,Xvi,Xvf,C,L);
            %  > 1.1.2. ---------------------------------------------------
            subplot(2,1,2);
            Fig_2_1D.SubPlot_2(msh,pde,Xvi,Xvf,C,L);           
        end
        %  > 1.1.1. -------------------------------------------------------
        function [] = SubPlot_1(msh,pde,Xvi,Xvf,C,L)
            hold on;
            P{1} = plot(msh.f.Xv,pde.e.f.f.a(:,1)               ,'-o','Color',C(1,:),'LineWidth',1.5,'MarkerFaceColor',C(1,:),'MarkerSize',3.5);
            P{2} = plot(msh.f.Xv,pde.e.f.f.a(:,2)               ,'-o','Color',C(2,:),'LineWidth',1.5,'MarkerFaceColor',C(2,:),'MarkerSize',3.5);
            P{3} = plot(msh.c.Xc,pde.e.c.c.a(:,1)               ,'-o','Color',C(3,:),'LineWidth',1.5,'MarkerFaceColor',C(3,:),'MarkerSize',3.5);
            P{4} = line([Xvi,Xvf],[pde.e.f.n.a(1,1),pde.e.f.n.a(1,1)],'Color',C(1,:),'Linewidth',1.0,'Linestyle','--');
            P{5} = line([Xvi,Xvf],[pde.e.f.n.a(1,2),pde.e.f.n.a(1,2)],'Color',C(2,:),'Linewidth',1.0,'Linestyle','-.');
            P{6} = line([Xvi,Xvf],[pde.e.c.n.a(1,1),pde.e.c.n.a(1,1)],'Color',C(3,:),'Linewidth',1.0,'Linestyle',':');
            %  > Legend.
            legend([P{1},P{2},P{3},P{4},P{5},P{6}],[L{1},L{2},L{3},L{4},L{5},L{6}],...
                'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',2);            
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{7},20,12);
        end
        %  > 1.1.2. -------------------------------------------------------
        function [] = SubPlot_2(msh,pde,Xvi,Xvf,C,L)
            hold on;
            P{1} = plot(msh.f.Xv,pde.e.f.f.f(:,1)               ,'-o','Color',C(1,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            P{2} = plot(msh.f.Xv,pde.e.f.f.f(:,2)               ,'-o','Color',C(2,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            P{3} = plot(msh.c.Xc,pde.e.c.c.c(:,1)               ,'-s','Color',C(3,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            P{4} = line([Xvi,Xvf],[pde.e.f.n.f(1,1),pde.e.f.n.f(1,1)],'Color',C(1,:),'Linewidth',1.0,'Linestyle','--');
            P{5} = line([Xvi,Xvf],[pde.e.f.n.f(1,2),pde.e.f.n.f(1,2)],'Color',C(2,:),'Linewidth',1.0,'Linestyle','-.');
            P{6} = line([Xvi,Xvf],[pde.e.c.n.c(1,1),pde.e.c.n.c(1,1)],'Color',C(3,:),'Linewidth',1.0,'Linestyle',':');
            %  > Legend.
            legend([P{1},P{2},P{3},P{4},P{5},P{6}],[L{8},L{9},L{10},L{11},L{12},L{13}],...
                'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',2);            
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{14},20,12);
        end
        % >> 1.2. ---------------------------------------------------------
        function [L] = Set_Labels()
            L{1}  = "$\tau_{f}^{\phantom{\nabla}\phi}$";
            L{2}  = "$\tau_{f}^{\nabla\phi}$";
            L{3}  = "$\tau_{c}^{\phantom{\nabla}\phi}$";
            L{4}  = "$|\!|\tau_{f}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{5}  = "$|\!|\tau_{f}^{\nabla\phi}|\!|_{1}$";
            L{6}  = "$|\!|\tau_{c}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{7}  = "$\textrm{Error magnitude}\left(\tau^{\phi}\right)$";
            L{8}  = "$e_{f}^{\phantom{\nabla}\phi}$";
            L{9}  = "$e_{f}^{\nabla\phi}$";
            L{10} = "$e_{c}^{\phantom{\nabla}\phi}$";
            L{11} = "$|\!|e_{f}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{12} = "$|\!|e_{f}^{\nabla\phi}|\!|_{1}$";
            L{13} = "$|\!|e_{c}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{14} = "$\textrm{Error magnitude}\left(e^{\phi}\right)$";
        end
    end
end      