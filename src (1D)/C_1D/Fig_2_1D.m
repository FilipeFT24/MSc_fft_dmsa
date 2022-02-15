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
            Xvi   = msh.f.Xv(1);
            Xvf   = msh.f.Xv(end);
            C1    = linspecer(9 ,'qualitative');
            C2    = linspecer(20,'sequential');
            L     = Fig_2_1D.Set_Labels();
            [m,n] = size(msh.s.xf);
            for i = 1:m
                for j = 1:n
                    p(j,i) = length(msh.s.xf{i,j});
                end
            end
            %  > 1.1.1. ---------------------------------------------------
            subplot(2,1,1);
            Fig_2_1D.SubPlot_1(msh,pde,Xvi,Xvf,C1,L);
            %  > 1.1.2. ---------------------------------------------------
            subplot(2,1,2);
            Fig_2_1D.SubPlot_2(msh,C2,p,m,n);           
        end
        %  > 1.1.1. -------------------------------------------------------
        function [] = SubPlot_1(msh,pde,Xvi,Xvf,C1,L)
            %  > Plot...
            hold on;
            P{1}  = plot(msh.f.Xv,pde.e.t.f(:,1)                ,'--o','Color',C1(1,:),'LineWidth',1.5,'MarkerFaceColor',C1(1,:),'MarkerSize',3.5);
            P{2}  = plot(msh.f.Xv,pde.e.t.f(:,2)                ,'-.s','Color',C1(2,:),'LineWidth',1.5,'MarkerFaceColor',C1(2,:),'MarkerSize',3.5);
            P{3}  = plot(msh.c.Xc,pde.e.t.c(:,1)                 ,':^','Color',C1(3,:),'LineWidth',1.5,'MarkerFaceColor',C1(3,:),'MarkerSize',3.5);
            P{4}  = line([Xvi,Xvf],[pde.e.t.n.f(1,1),pde.e.t.n.f(1,1)],'Color',C1(1,:),'Linewidth',0.5,'Linestyle','--');
            P{5}  = line([Xvi,Xvf],[pde.e.t.n.f(1,2),pde.e.t.n.f(1,2)],'Color',C1(2,:),'Linewidth',0.5,'Linestyle','-.');
            P{6}  = line([Xvi,Xvf],[pde.e.t.n.c(1,1),pde.e.t.n.c(1,1)],'Color',C1(3,:),'Linewidth',0.5,'Linestyle',':' );
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',2);
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{7},20,12);
        end
        %  > 1.1.2. -------------------------------------------------------
        function [] = SubPlot_2(msh,C2,p,m,n)
            %  > Auxiliary variables.
            N = 5;
            Y = [2/N,1/N,1];  
            M = ['o','s'];
            %  > Plot...
            hold on;
            for i = 1:m
                plot(msh.f.Xv,repelem(Y(i),msh.f.NF),'-','Color','k','MarkerFaceColor','k','MarkerSize',10.0,'Linewidth',0.5);
                pm{1} = unique(p(:,i));
            end
            hold on;
            for i = 1:m+1
                if i ~= m+1
                    for j = 1:n
                        plot(msh.f.Xv(j),Y(i),M(i),'Color',C2(p(j,i),:),'MarkerFaceColor',C2(p(j,i),:),'MarkerSize',4.5);
                    end
                else
                    plot(0.5,(N-1)/N,'o','Color','w','MarkerFaceColor','w','MarkerSize',4.5);
                end 
            end
            %  > Legend.
            for i = 1:m
                pm{i} = unique(p(:,i));
            end
            i = 1;
            for j = 1:length(pm{i})
                Pj  {j} = plot(NaN,NaN,M(i),'Color',C2(pm{i}(j),:),'MarkerFaceColor',C2(pm{i}(j),:),'MarkerSize',3.5);
                Lj{1,j} = join(['$',num2str(pm{i}(j)),'\phantom{:}$']);
            end
            i = 2;
            for k = 1:length(pm{i})
                Pk  {k} = plot(NaN,NaN,M(i),'Color',C2(pm{i}(k),:),'MarkerFaceColor',C2(pm{i}(k),:),'MarkerSize',3.5);
                Lk{1,k} = join(['$',num2str(pm{i}(k)),'\phantom{:}$']);
            end
            legendflex([Pk{:}],Lk,'buffer',[-0.075,-0.025],'bufferunit', 'normalized',...
                'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',1);
            legendflex([Pj{:}],Lj,'buffer',[-0.010,-0.025],'bufferunit', 'normalized',...
                'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',1);
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(false,true,false,msh.f.Xv,10,"$x$","$y$",20,12);
            ax = gca; box on; ylim([0,1]); set(gca,'YTickLabel',[]); yticks([0,1]);
        end
        % >> 1.2. ---------------------------------------------------------
        function [L] = Set_Labels()
            L{1}  = "$|\tau_{f}^{\phantom{\nabla}\phi}|$";
            L{2}  = "$|\tau_{f}^{\nabla\phi}|$";
            L{3}  = "$|\tau_{c}^{\phantom{\nabla}\phi}|$";
            L{4}  = "$|\!|\tau_{f}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{5}  = "$|\!|\tau_{f}^{\nabla\phi}|\!|_{1}$";
            L{6}  = "$|\!|\tau_{c}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{7}  = "$\textrm{Error magnitude}\left(\tau^{\phi}\right)$";
            L{8}  = "$|e_{f}^{\phantom{\nabla}\phi}|$";
            L{9}  = "$|e_{f}^{\nabla\phi}|$";
            L{10} = "$|e_{c}^{\phantom{\nabla}\phi}|$";
            L{11} = "$|\!|e_{f}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{12} = "$|\!|e_{f}^{\nabla\phi}|\!|_{1}$";
            L{13} = "$|\!|e_{c}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{14} = "$\textrm{Error magnitude}\left(e^{\phi}\right)$";
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = SubPlot_3(msh,pde,Xvi,Xvf,C2,L)
            hold on;
            P{1} = plot(msh.f.Xv,pde.e.f.f(:,1)             ,'-o','Color',C2(1,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            P{2} = plot(msh.f.Xv,pde.e.f.f(:,2)             ,'-o','Color',C2(2,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            P{3} = plot(msh.c.Xc,pde.e.c.c(:,1)             ,'-s','Color',C2(3,:),'LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',3.5);
            P{4} = line([Xvi,Xvf],[pde.e.f.n(1,1),pde.e.f.n(1,1)],'Color',C2(1,:),'Linewidth',1.0,'Linestyle','--');
            P{5} = line([Xvi,Xvf],[pde.e.f.n(1,2),pde.e.f.n(1,2)],'Color',C2(2,:),'Linewidth',1.0,'Linestyle','-.');
            P{6} = line([Xvi,Xvf],[pde.e.c.n(1,1),pde.e.c.n(1,1)],'Color',C2(3,:),'Linewidth',1.0,'Linestyle',':' );
            %  > Legend.
            legend([P{1},P{2},P{3},P{4},P{5},P{6}],[L{8},L{9},L{10},L{11},L{12},L{13}],...
                'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',2);            
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{14},20,12);
        end
    end
end      