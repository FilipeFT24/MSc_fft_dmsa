classdef Fig_4_1D
    methods (Static)         
        %% > Wrap-up Fig_4 (1D).
        function [] = WrapUp_Fig_4_1D(Plot_4,Exp_4,Fig,Xv,n_LH,XN,E1,E2)
            if Plot_4               
                %  > Figure 1.
                %figure(Fig(1)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                %Fig_4_1D.Plot_1(Xv,XN,E1,E2,n_LH);
                %  > Figure 2.
                figure(Fig(2)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_4_1D.Plot_2(Xv,XN,E1,E2);
                %  > Export as .pdf.
                if Exp_4
                    Fig_Tools_1D.Export_PDF('Fig_4','../[Figures]/[1D]/Fig_4');
                end
            end
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [] = Plot_1(Xv,E1,E2,n_LH)
            %  > 1.1.0. ---------------------------------------------------
            C  = linspecer(9,'qualitative');
            L  = Fig_4_1D.Set_Labels_1();
            %  > 1.1.1. ---------------------------------------------------
            subplot(1,2,1);
            hold on;
            for i = 1:size(E2,2)
                P1(i) = plot(Xv,E1{i+1}.f.f(:,1)./E2{i}.f.f(:,1),'-s','Color',C(i,:),'LineWidth',1.25,'MarkerFaceColor','w','MarkerSize',3.0);
                L1(i) = convertCharsToStrings(num2str(n_LH(i+1)));
            end
            legend(P1,L1,'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',1);
            Fig_Tools_1D.ChangeLook_1D(true,true,Xv,10,"$x$",L{1},20,12);
            %  > 1.1.2. ---------------------------------------------------
            subplot(1,2,2);
            hold on;
            for i = 1:size(E2,2)
                P2(i) = plot(Xv,E1{i+1}.f.f(:,2)./E2{i}.f.f(:,2),'-^','Color',C(i,:),'LineWidth',1.25,'MarkerFaceColor','w','MarkerSize',3.0);
                L2(i) = convertCharsToStrings(num2str(n_LH(i+1)));
            end
            legend(P2,L2,'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',1);
            Fig_Tools_1D.ChangeLook_1D(true,true,Xv,10,"$x$",L{2},20,12);
        end
        % >> 1.2. ---------------------------------------------------------
        function [L] = Set_Labels_1()
            L{1} = "$\zeta_{f}^{\phi}$";
            L{2} = "$\zeta_{f}^{\nabla\phi}$";
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_2(Xv,XN,E1,E2)
            %  > 1.1.0. ---------------------------------------------------
            C    = linspecer(9,'qualitative');
            L    = Fig_4_1D.Set_Labels_2();
            %  > 1.1.1. ---------------------------------------------------
            subplot(1,2,1);
            hold on;
            P{1} = plot(Xv,E2{1}.f.f(:,1)./max(E2{1}.f.f(:,1)),'-s','Color',C(1,:),'LineWidth',1.25,'MarkerFaceColor','w','MarkerSize',3.0);
            P{2} = plot(Xv,E1{2}.f.f(:,1)./max(E1{2}.f.f(:,1)),'-s','Color',C(2,:),'LineWidth',1.25,'MarkerFaceColor','w','MarkerSize',3.0);
            Fig_Tools_1D.ChangeLook_1D(true,true,Xv,10,"$x$",L{5},20,12);
            legend([P{1},P{2}],[L{1},L{2}],'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',1);
            %  > 1.1.2. ---------------------------------------------------
            subplot(1,2,2);
            hold on;
            P{3} = plot(Xv,E2{1}.f.f(:,2)./max(E2{1}.f.f(:,2)),'-^','Color',C(1,:),'LineWidth',1.25,'MarkerFaceColor','w','MarkerSize',3.0);
            P{4} = plot(Xv,E1{2}.f.f(:,2)./max(E1{2}.f.f(:,2)),'-^','Color',C(2,:),'LineWidth',1.25,'MarkerFaceColor','w','MarkerSize',3.0);
            legend([P{3},P{4}],[L{3},L{4}],'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',1);
            Fig_Tools_1D.ChangeLook_1D(true,true,Xv,10,"$x$",L{5},20,12);
        end
        % >> 1.2. ---------------------------------------------------------
        function [L] = Set_Labels_2()
            L{1} = "$E_{f}^{\phi}$";
            L{2} = "$\epsilon_{f}^{\phi}$";
            L{3} = "$E_{f}^{\nabla\phi}$";
            L{4} = "$\epsilon_{f}^{\nabla\phi}$";
            L{5} = "$\textrm{Error magnitude}$";
        end
    end
end      