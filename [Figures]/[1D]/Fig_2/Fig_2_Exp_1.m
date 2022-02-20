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
                %if Exp_2
                %    Fig_Tools_1D.Export_PDF('ET_4_1_1','../[Figures]/[1D]/Fig_2');
                %end
                figure(Fig(2)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_2(msh,pde);
                %  > Export as .pdf.
                if Exp_2
                    Fig_Tools_1D.Export_PDF('ET_4_2 (3)','../[Figures]/[1D]/Fig_2');
                end
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
            %subplot(2,1,1);
            Fig_2_1D.SubPlot_1(msh,pde,Xvi,Xvf,C1,L);
            %  > 1.1.2. ---------------------------------------------------
            %subplot(2,1,2);
            %Fig_2_1D.SubPlot_3(msh,C2,p,m,n);           
        end
        %  > 1.1.1. -------------------------------------------------------
        function [] = SubPlot_1(msh,pde,Xvi,Xvf,C1,L)
            %  > Plot...
            m = size(pde.e.t.f,2)+1;
            hold on;
            for j = 2:m
                if j == 2%m
                    P{1}   = plot(msh.f.Xv,pde.e.t.f(:,j),':o','Color',C1(j,:),'LineWidth',3.0,'MarkerFaceColor',C1(j,:),'MarkerSize',5.0);
                    P{2} = line([Xvi,Xvf],[pde.e.t.n.f(1,j),pde.e.t.n.f(1,j)],'Color',C1(j,:),'Linewidth',2.5,'Linestyle','-.');
                else
                    %P{j+m} = line([Xvi,Xvf],[pde.e.t.n.t(1,1),pde.e.t.n.t(1,1)],'Color',C1(j,:),'Linewidth',1.0,'Linestyle','--');
                end
            end
            P{3} = plot(msh.c.Xc,abs(pde.e.t.c(:,1)),':o','Color',C1(j,:),'LineWidth',3.0,'MarkerFaceColor',C1(j,:),'MarkerSize',5.0);
            P{4} = line([Xvi,Xvf],[pde.e.t.n.c(1,1),pde.e.t.n.c(1,1)],'Color',C1(j,:),'Linewidth',2.5,'Linestyle','-.');
            
            
            
            
            
            legend([P{:}],[L{2},L{4},L{6},L{7},L{8}],...
                'Interpreter','latex','Location','Northeast','FontSize',25,'NumColumns',2);
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{8},30,30);
            ylim([10^-10,10^0]);
            set(gca,'YScale','log');
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
            %subplot(2,1,1);
            Fig_2_1D.SubPlot_2(msh,pde,Xvi,Xvf,C1,L);
            %  > 2.1.2. ---------------------------------------------------
            %subplot(2,1,2);
            %Fig_2_1D.SubPlot_3(msh,C2,p,m,n);           
        end
        %  > 2.1.1. -------------------------------------------------------
        function [] = SubPlot_2(msh,pde,Xvi,Xvf,C1,L)
            hold on;
            P{1} = plot(msh.c.Xc,pde.e.c.c(:,1)             ,'-s','Color',C1(4,:),'LineWidth',3.0,'MarkerFaceColor','w','MarkerSize',5.0);
            P{2} = line([Xvi,Xvf],[pde.e.c.n(1,1),pde.e.c.n(1,1)],'Color',C1(4,:),'Linewidth',2.5,'Linestyle',':' );
            
            D1 = [2.00804903165388e-12;3.39590099137216e-12;6.82667385015941e-12;1.51263977120562e-11;3.47760531783744e-11;8.03002661216016e-11;1.83503469680452e-10;4.12415368187323e-10;9.09150256052256e-10;1.96357192972048e-09;4.15276946811184e-09;8.59784344410093e-09;1.74232167180206e-08;3.45539999867730e-08;6.70577731856728e-08;1.27329871149991e-07;2.36531240172670e-07;4.29799402758337e-07;7.63832865241820e-07;1.32744070313381e-06;2.25546849571417e-06;3.74605266826770e-06;6.08030426262789e-06;9.64219361414222e-06;1.49345925179584e-05;2.25852709480265e-05;3.33344979047091e-05;4.79943758475545e-05;6.73700010399500e-05;9.21349676730760e-05;0.000122659488567888;0.000158798895283939;0.000199663050711126;0.000243401606350349;0.000287053100687093;0.000326513608421355;0.000356678675877492;0.000371797123458362;0.000366045702174911;0.000366203366004625;0.000366324008211794;0.000366105601435018;0.000365499401928648;0.000364515635608242;0.000363222385522777;0.000361742434581980;0.000360239790724415;0.000358897045970474;0.000357887222238174;0.000357345333931058;0.000357345333930836;0.000357887222237396;0.000358897045969364;0.000360239790723083;0.000361742434580536;0.000363222385521889;0.000364515635607798;0.000365499401927871;0.000366105601433298;0.000366324008210017;0.000366203366003071;0.000366045702173412;0.000371797123456613;0.000356678675875494;0.000326513608419315;0.000287053100684817;0.000243401606348129;0.000199663050708870;0.000158798895281656;0.000122659488565754;9.21349676710412e-05;6.73700010380210e-05;4.79943758456957e-05;3.33344979029379e-05;2.25852709463265e-05;1.49345925163252e-05;9.64219361257674e-06;6.08030426113018e-06;3.74605266683579e-06;2.25546849434732e-06;1.32744070183311e-06;7.63832864007259e-07;4.29799401590488e-07;2.36531239071444e-07;1.27329870115575e-07;6.70577722180465e-08;3.45539990858866e-08;1.74232158838685e-08;8.59784267668183e-09;4.15276876742392e-09;1.96357129576480e-09;9.09149688828656e-10;4.12414867695858e-10;1.83503035921167e-10;8.02998990945381e-11;3.47757528834925e-11;1.51261641493688e-11;6.82650701967002e-12;3.39580089307889e-12;2.00801566555624e-12];
            DN = [0.000134473830667821;0.000134473830667821;0.000371797123458362];
            
            
            P{3} = plot(msh.c.Xc,D1(:,1)           ,'-s','Color',C1(1,:),'LineWidth',3.0,'MarkerFaceColor','w','MarkerSize',5.0);
            P{4} = line([Xvi,Xvf],[DN(1,1),DN(1,1)],'Color',C1(1,:),'Linewidth',2.5,'Linestyle',':' );
            
            %  > Legend.
            %legend([P{:}],[L{:}],...
            %    'Interpreter','latex','Location','Northeast','FontSize',25,'NumColumns',2);
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{5},30,30);
            
            
            Pk  {1} = plot(NaN,NaN,'-s','Color','k','MarkerFaceColor','w','MarkerSize',5.0);
            Pk  {2} = plot(NaN,NaN     ,'Color','k','Linewidth',2.5,'Linestyle',':');
            Pj  {1} = plot(NaN,NaN,'-s','Color',C1(1,:),'MarkerFaceColor',C1(1,:),'MarkerSize',5.0);
            Pj  {2} = plot(NaN,NaN,'-s','Color',C1(4,:),'MarkerFaceColor',C1(4,:),'MarkerSize',5.0);
            
            Lk{1,1} = '$|e_{c}^{\phantom{\nabla}\phi}|$';
            Lk{1,2} = '$|\!|e_{c}^{\phantom{\nabla}\phi}|\!|_{1}$';
            Lj{1,1} = '$_{\mathcal{D}_{1}}$';
            Lj{1,2} = '$_{\mathcal{D}_{1}+\mathcal{D}_{2}}$';
            
            
            
            legendflex([Pj{:}],Lj,'buffer',[-0.010,-0.025],'bufferunit', 'normalized',...
                'Interpreter','latex','Location','Northeast','FontSize',25,'NumColumns',2);
            legendflex([Pk{:}],Lk,'buffer',[-0.130,-0.025],'bufferunit', 'normalized',...
                'Interpreter','latex','Location','Northeast','FontSize',25,'NumColumns',2);
            
        end
        % >> 2.2. ---------------------------------------------------------
        function [L] = Set_Labels_2()
            L{1} = "$|e_{c}^{\phantom{\nabla}\phi}|$";
            L{2} = "$|\!|e_{c}^{\phantom{\nabla}\phi}|\!|_{1}$";
            
            L{3} = "$\mathcal{D}_{1}$";
            L{4} = "$\mathcal{D}_{1}+\mathcal{D}_{2}$";
            L{5} = "$\textrm{Error magnitude}, |e^{\phi}|$";
        end
    end
end      