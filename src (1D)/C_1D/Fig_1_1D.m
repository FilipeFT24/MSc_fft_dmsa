classdef Fig_1_1D
    methods (Static)
        %% > #1.
        function [] = WrapUp_Fig_1_1_1D(flag_1,msh,pde)
            if flag_1
                %  > Plot/Export.
                Exp = false;
                F_1 = "FN_1";
                F_2 = "FN_2";
                F_3 = "FN_3";
                D_1 = "../[Figures]/[1D]/Fig_2";
                D_2 = "../[Figures]/[1D]/Fig_2";
                D_3 = "../[Figures]/[1D]/Fig_2";
                %  > Properties.
                Fig = [1,2,3];
                fig = Fig_1_1D.Set_fig_1(Exp);
                l_1 = true;
                l_2 = true;
                l_3 = true;
                
                if ~Exp
                    %  > 1/2.
                    figure(Fig(1)); set(gcf,'Units','pixels','Position',[150,100,1250,600]);
                    subplot(1,2,1);
                    Fig_1_1D.Plot_1_1(fig,l_1,msh,pde);
                    subplot(1,2,2);
                    Fig_1_1D.Plot_1_2(fig,l_2,msh,pde);
                    %  > 3.
                    %  figure(Fig(3)); set(gcf,'Units','pixels','Position',[500,100,600,600]);
                    %  Fig_1_1D.Plot_1_3(fig,msh);
                else
                    %  > 1.
                    figure(Fig(1)); set(gcf,'Units','pixels','Position',[350,100,850,600]);
                    Fig_1_1D.Plot_1_1(fig,l_1,msh,pde);
                    Fig_Tools_1D.Export_PDF(F_1,D_1);
                    %  > 2.
                    figure(Fig(2)); set(gcf,'Units','pixels','Position',[350,100,850,600]);
                    Fig_1_1D.Plot_1_2(fig,l_2,msh,pde);
                    Fig_Tools_1D.Export_PDF(F_2,D_2);
                end
            end
        end
        % >> 1. -----------------------------------------------------------
        function [] = Plot_1_1(fig,l_1,msh,pde)
            %  > Auxiliary variables.
            C    = linspecer(9,'qualitative');
            L    = Fig_1_1D.Set_Labels_1_1();
            m    = size(pde.e.t.f,2);
            Xv   = [msh.f.Xv(1),msh.f.Xv(msh.f.NF)];
            trsh = 10e-12;
            
            hold on;
            for j = 1:m
                %  > tau_f.
                i        = pde.e.t.f_abs(:,j) > trsh;
                P{j}     = plot(msh.f.Xv,pde.e.t.f_abs(:,j)              ,':o','Color',C(j,:),'LineWidth',fig.LW_1,'MarkerFaceColor',C(j,:),'MarkerSize',fig.MS_1);
                P{j+m+1} = line(Xv,[pde.e.t.n_abs.f(1,j),pde.e.t.n_abs.f(1,j)],'Color',C(j,:),'Linewidth',fig.LW_2,'Linestyle','-');
                if l_1
                    y_min(j) = min(pde.e.t.f_abs(i,j));
                    y_max(j) = max(pde.e.t.f_abs(i,j));
                end
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',fig.FT_1,'NumColumns',2);
            %  > Axis.
            if l_1
                set(gca,'YScale','log');
                ylim([10.^(ceil(log10(min(y_min)))-1),...
                      10.^(ceil(log10(max(y_max)))+1)]);
            end
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{end},fig.FT_2,fig.FT_3);
        end
        % >> 2. -----------------------------------------------------------
        function [] = Plot_1_2(fig,l_2,msh,pde)
            %  > Auxiliary variables.
            C    = linspecer(9,'qualitative');
            L    = Fig_1_1D.Set_Labels_1_2();
            Xv   = [msh.f.Xv(1),msh.f.Xv(msh.f.NF)];
            n    = 2;
            
            hold on;
            for i = 1:n
                if i ~= n
                    %  > e_c.
                    P{i}     = plot(msh.c.Xc,pde.e.c.c_abs(:,1)              ,':o','Color',C(i,:),'LineWidth',fig.LW_1,'MarkerFaceColor',C(i,:),'MarkerSize',fig.MS_1);
                    P{i+n}   = line(Xv,[pde.e.c.n_abs(1,1),pde.e.c.n_abs(1,1)]    ,'Color',C(i,:),'Linewidth',fig.LW_2,'Linestyle','-');
                    y_min(i) = min(pde.e.c.c_abs(:,1));
                    y_max(i) = max(pde.e.c.c_abs(:,1));
                else
                    %  > tau_c.
                    P{i}     = plot(msh.c.Xc,pde.e.t.c_abs(:,1)            , '-.o','Color',C(i,:),'LineWidth',fig.LW_1,'MarkerFaceColor',C(i,:),'MarkerSize',fig.MS_1);
                    P{i+n}   = line(Xv,[pde.e.t.n_abs.c(1,1),pde.e.t.n_abs.c(1,1)],'Color',C(i,:),'Linewidth',fig.LW_2,'Linestyle','-');
                    y_min(i) = min(pde.e.t.c_abs(:,1));
                    y_max(i) = max(pde.e.t.c_abs(:,1));
                    
                end
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',fig.FT_1,'NumColumns',2);
            %  > Axis.
            if l_2
                set(gca,'YScale','log');
                ylim([10.^(ceil(log10(min(y_min)))-1),...
                    10.^(ceil(log10(max(y_max)))+1)]);
            end
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{end},fig.FT_2,fig.FT_3);
        end
        % >> 3. -----------------------------------------------------------
        function [] = Plot_1_3(fig,msh)
            %  > Auxiliary variables.
            C = linspecer(9,'qualitative');
            L = Fig_1_1D.Set_Labels_1_3();
            m = msh.f.NF;
            n = size(msh.s.stl.p,2);
            M = ["-o",":d"];
            
            hold on;
            for i = 1:n
                for j = 1:m
                    p(j,i) = A_2_1D.Compute_p(msh.s.stl.p{i}(j),msh.s.stl.t{i}(j));
                end
                P{i} = stem(msh.f.Xv,p(:,i),M(i),'Color',C(i,:),'Linewidth',fig.LW_2);
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',fig.FT_1,'NumColumns',2);
            %  > Axis.
            ax = get(gca,'yTick'); yticks(unique(round(ax)));
            ylim([0,max(p,[],'all')+0.5]);
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{end},fig.FT_2,fig.FT_3);
        end
        % >> 4. -----------------------------------------------------------
        %  > 4.1. ---------------------------------------------------------
        function [L] = Set_Labels_1_1()
            L{1} = "$|\tau_{f}^{\phantom{\nabla}\phi}|$";
            L{2} = "$|\tau_{f}^{\nabla\phi}|$";
            L{3} = "$|\tau_{f}^{\phantom{\nabla\phi}}|$";
            L{4} = "$|\!|\tau_{f}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{5} = "$|\!|\tau_{f}^{\nabla\phi}|\!|_{1}$";
            L{6} = "$|\!|\tau_{f}^{\phantom{\nabla\phi}}|\!|_{1}$";
            L{7} = "$\textrm{Error magnitude}$";
        end
        %  > 4.2. ---------------------------------------------------------
        function [L] = Set_Labels_1_2()
            L{1} = "$\phantom{|\!}|e_{c}|$";
            L{2} = "$\phantom{|\!}|\tau_{c}|$";
            L{3} = "$|\!|e_{c}|\!|_{1}$";
            L{4} = "$|\!|\tau_{c}|\!|_{1}$";
            L{5} = "$\textrm{Error magnitude}$";
        end
        %  > 4.3. ---------------------------------------------------------
        function [L] = Set_Labels_1_3()
            L{1} = "$\phi_{f}$";
            L{2} = "$\nabla\phi_{f}$";
            L{3} = "$\textrm{Method's order}$";
        end
        %  > 4.4. ---------------------------------------------------------
        function [fig] = Set_fig_1(Exp)
            if ~Exp
                fig.LW_1 =  2.0;
                fig.LW_2 =  1.5;
                fig.MS_1 =  3.0;
                fig.FT_1 = 11.5;
                fig.FT_2 = 12.5;
                fig.FT_3 = 12.5;
            else
                fig.LW_1 =  3.0;
                fig.LW_2 =  2.5;
                fig.MS_1 =  5.0;
                fig.FT_1 = 20.0;
                fig.FT_2 = 30.0;
                fig.FT_3 = 30.0;
            end
        end
        
        %% > #2.
        function [] = WrapUp_Fig_1_2_1D(flag_2,msh,pde)
            if flag_2
                %  > Properties.
                Fig = [4,5];
                fig = Fig_1_1D.Set_fig_2();
                
                figure(Fig(1)); set(gcf,'Units','pixels','Position',[150,100,1250,600]);
                subplot(1,2,1);
                Fig_1_1D.Plot_2_1(fig,msh,pde);
                subplot(1,2,2);
                Fig_1_1D.Plot_2_2(fig,msh,pde);
            end
        end
        % >> 1. -----------------------------------------------------------
        function [] = Plot_2_1(fig,msh,pde)
            %  > Auxiliary variables.
            C  = linspecer(9,'qualitative');
            L  = Fig_1_1D.Set_Labels_2_1();
            m  = size(pde.e.t.f,2);
            Xv = [msh.f.Xv(1),msh.f.Xv(msh.f.NF)];
            
            %  > tau_f.
            hold on;
            for j = 1:m
                P{j}     = plot(msh.f.Xv,pde.e.t.f(:,j)          ,':o','Color',C(j,:),'LineWidth',fig.LW_1,'MarkerFaceColor',C(j,:),'MarkerSize',fig.MS_1);
                P{j+m+1} = line(Xv,[pde.e.t.n.f(1,j),pde.e.t.n.f(1,j)],'Color',C(j,:),'Linewidth',fig.LW_2,'Linestyle','-');
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',fig.FT_1,'NumColumns',2);
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{end},fig.FT_2,fig.FT_3);
        end
        % >> 2. -----------------------------------------------------------
        function [] = Plot_2_2(fig,msh,pde)
            %  > Auxiliary variables.
            C  = linspecer(9,'qualitative');
            L  = Fig_1_1D.Set_Labels_2_2();
            Xv = [msh.f.Xv(1),msh.f.Xv(msh.f.NF)];
            n  = 2;
            
            hold on;
            for i = 1:n
                if i ~= n
                    %  > e_c.
                    P{i}   = plot(msh.c.Xc,pde.e.c.c(:,1)          ,':o','Color',C(i,:),'LineWidth',fig.LW_1,'MarkerFaceColor',C(i,:),'MarkerSize',fig.MS_1);
                    P{i+n} = line(Xv,[pde.e.c.n(1,1),pde.e.c.n(1,1)]    ,'Color',C(i,:),'Linewidth',fig.LW_2,'Linestyle','-');
                else
                    %  > tau_c.
                    P{i}   = plot(msh.c.Xc,pde.e.t.c(:,1)          ,':o','Color',C(i,:),'LineWidth',fig.LW_1,'MarkerFaceColor',C(i,:),'MarkerSize',fig.MS_1);
                    P{i+n} = line(Xv,[pde.e.t.n.c(1,1),pde.e.t.n.c(1,1)],'Color',C(i,:),'Linewidth',fig.LW_2,'Linestyle','-');
                end
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',fig.FT_1,'NumColumns',2);
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",L{end},fig.FT_2,fig.FT_3);
        end
        % >> 3. -----------------------------------------------------------
        %  > 3.1. ---------------------------------------------------------
        function [L] = Set_Labels_2_1()
            L{1} = "$\tau_{f}^{\phantom{\nabla}\phi}$";
            L{2} = "$\tau_{f}^{\nabla\phi}$";
            L{3} = "$\tau_{f}^{\phantom{\nabla\phi}}$";
            L{4} = "$|\!|\tau_{f}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{5} = "$|\!|\tau_{f}^{\nabla\phi}|\!|_{1}$";
            L{6} = "$|\!|\tau_{f}^{\phantom{\nabla\phi}}|\!|_{1}$";
            L{7} = "$\textrm{Error magnitude}$";
        end
        %  > 3.2. ---------------------------------------------------------
        function [L] = Set_Labels_2_2()
            L{1} = "$\phantom{|\!|}e_{c}$";
            L{2} = "$\phantom{|\!|}\tau_{c}$";
            L{3} = "$|\!|e_{c}|\!|_{1}$";
            L{4} = "$|\!|\tau_{c}|\!|_{1}$";
            L{5} = "$\textrm{Error magnitude}$";
        end
        %  > 3.3. ---------------------------------------------------------
        function [fig] = Set_fig_2()
            fig.LW_1 =  2.0;
            fig.LW_2 =  1.5;
            fig.MS_1 =  3.0;
            fig.FT_1 = 11.5;
            fig.FT_2 = 12.5;
            fig.FT_3 = 12.5;
        end
    end
end