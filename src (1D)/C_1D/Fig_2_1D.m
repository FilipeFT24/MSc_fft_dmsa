classdef Fig_2_1D
    methods (Static)         
        %% > Wrap-up Fig_2 (1D).
        function [] = WrapUp_Fig_2_1D(Plot_2,Exp_2,Fig,msh,pde)
            if Plot_2               
                %  > Figure 1.
                figure(Fig(1)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_1(msh,pde);
                %  > Figure 2.
                %  > Export as .pdf.
                if Exp_2
                    Fig_Tools_1D.Export_PDF('ET_4_1 (3)','../[Figures]/[1D]/Fig_2');
                end
                %  figure(Fig(2)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                %  Fig_2_1D.Plot_2(msh,pde);
                %  %  > Export as .pdf.
                %  if Exp_2
                %      Fig_Tools_1D.Export_PDF('ET_4_2 (3)','../[Figures]/[1D]/Fig_2');
                %  end
            end
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [] = Plot_1(msh,pde)
            %  > Auxiliary variables.
            Xvi   = msh.f.Xv(1);
            Xvf   = msh.f.Xv(end);
            C1    = linspecer(9 ,'qualitative');
            C2    = linspecer(20,'sequential');
            L     = Fig_2_1D.Set_Labels_1();
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
            Fig_2_1D.SubPlot_3(msh,C2,p,m,n);           
        end
        %  > 1.1.1. -------------------------------------------------------
        function [] = SubPlot_1(msh,pde,Xvi,Xvf,C1,L)
            %  > Plot...
            m = size(pde.e.t.f,2);
            hold on;
            for j = 1:m
                P{j}   = plot(msh.f.Xv,pde.e.t.f(:,j)                 ,':o','Color',C1(j,:),'LineWidth',2.0,'MarkerFaceColor',C1(j,:),'MarkerSize',3.5);
                P{j+m} = line([Xvi,Xvf],[pde.e.t.n.f(1,j),pde.e.t.n.f(1,j)],'Color',C1(j,:),'Linewidth',1.5,'Linestyle','-.');
            end
            
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',2);
            %  > Axis.
            set(gca,'YScale','log');
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{8},15,15); 
        end
        %  > 1.1.2. -------------------------------------------------------
        function [] = SubPlot_3(msh,C2,p,m,n)
            %  > Auxiliary variables.
            N = 9;
            Y = [2/N,1/N,1];  
            M = ["-o","o","-s","s"];
            %  > Plot...
            hold on;
            for i = 1:m
                plot(msh.f.Xv,repelem(Y(i),msh.f.NF),'-','Color','k','MarkerFaceColor','k','MarkerSize',10.0,'Linewidth',0.1);
                pm{1} = unique(p(:,i));
            end
            hold on;
            for i = 1:m+1
                if i ~= m+1
                    for j = 1:n
                        plot(msh.f.Xv(j),Y(i),M(2*i-1),'Color',C2(p(j,i),:),'MarkerFaceColor',C2(p(j,i),:),'MarkerSize',4.5);
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
            for j = 1:length(pm{i})+1
                if j == 1
                    Pj  {j} = plot(NaN,NaN,M(i),'Color','k','MarkerFaceColor','w','MarkerSize',3.5);
                    Lj{1,j} = '${\phi}_{f}\phantom{:}$';
                else
                    Pj  {j} = plot(NaN,NaN,M(i+1),'Color',C2(pm{i}(j-1),:),'MarkerFaceColor',C2(pm{i}(j-1),:),'MarkerSize',3.5);
                    Lj{1,j} = join(['$',num2str(pm{i}(j-1)),'\phantom{:}$']);
                end
            end
            i = 2;
            for k = 1:length(pm{i})+1
                if k == 1
                    Pk  {k} = plot(NaN,NaN,M(i+1),'Color','k','MarkerFaceColor','w','MarkerSize',3.5);
                    Lk{1,k} = '${\nabla\phi}_{f}\phantom{:}$';
                else
                    Pk  {k} = plot(NaN,NaN,M(i+2),'Color',C2(pm{i}(k-1),:),'MarkerFaceColor',C2(pm{i}(k-1),:),'MarkerSize',3.5);
                    Lk{1,k} = join(['$',num2str(pm{i}(k-1)),'\phantom{:}$']);
                end
            end
            legendflex([Pj{:}],Lj,'buffer',[-0.090,-0.025],'bufferunit', 'normalized',...
                'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',2);
            legendflex([Pk{:}],Lk,'buffer',[-0.010,-0.025],'bufferunit', 'normalized',...
                'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',2);
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(false,true,false,msh.f.Xv,10,"$x$","$y$",20,12);
            ax = gca; box on; ylim([0,1]); set(gca,'YTickLabel',[]); yticks([0,1]);
        end
        % >> 1.2. ---------------------------------------------------------
        function [L] = Set_Labels_1()
            L{1} = "$|\tau_{f}^{\phantom{\nabla}\phi}|$";
            L{2} = "$|\tau_{f}^{\nabla\phi}|$";
            L{3} = "$|\!|\tau_{f}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{4} = "$|\!|\tau_{f}^{\nabla\phi}|\!|_{1}$";
            L{5} = "$|\!|\bar\tau_{f}^{\phantom{\nabla\phi}}|\!|_{1}$";
            L{6} = "$|\tau_{c}|$";
            L{7} = "$|\!|\tau_{c}|\!|_{1}$";
            L{8} = "$\textrm{Error magnitude}, |\tau^{\phi}|$";
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_2(msh,pde)
            %  > Auxiliary variables.
            Xvi   = msh.f.Xv(1);
            Xvf   = msh.f.Xv(end);
            C1    = linspecer(9 ,'qualitative');
            C2    = linspecer(20,'sequential');
            L     = Fig_2_1D.Set_Labels_2();
            [m,n] = size(msh.s.xf);
            for i = 1:m
                for j = 1:n
                    p(j,i) = length(msh.s.xf{i,j});
                end
            end
            %  > 2.1.1. ---------------------------------------------------
            subplot(2,1,1);
            Fig_2_1D.SubPlot_2(msh,pde,Xvi,Xvf,C1,L);
            %  > 2.1.2. ---------------------------------------------------
            subplot(2,1,2);
            Fig_2_1D.SubPlot_3(msh,C2,p,m,n);           
        end
        %  > 2.1.1. -------------------------------------------------------
        function [] = SubPlot_2(msh,pde,Xvi,Xvf,C1,L)
            hold on;
            P{1} = plot(msh.c.Xc,pde.e.c.c(:,1)             ,'-s','Color',C1(1,:),'LineWidth',2.0,'MarkerFaceColor','w','MarkerSize',3.5);
            P{2} = line([Xvi,Xvf],[pde.e.c.n(1,1),pde.e.c.n(1,1)],'Color',C1(1,:),'Linewidth',1.5,'Linestyle',':' );
            %  > Legend.
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',2);
            %  > Axis.
            set(gca,'YScale','log');
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{3},15,15);
        end
        % >> 2.2. ---------------------------------------------------------
        function [L] = Set_Labels_2()
            L{1} = "$|e_{c}^{\phantom{\nabla}\phi}|$";
            L{2} = "$|\!|e_{c}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{3} = "$\textrm{Error magnitude}, |e^{\phi}|$";
        end
    end
end      