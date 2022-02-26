classdef Fig_2_1D
    methods (Static)
        %% > #1.
        function [] = WrapUp_Fig_2_1_1D(msh,pde,p,ttm)
            %  > Plot/Export.
            Exp = false;
            F_1 = "FN_1";
            F_2 = "FN_2";
            D_1 = "../[Figures]/[1D]/Fig_2";
            D_2 = "../[Figures]/[1D]/Fig_2";
            %  > Properties.
            fig = Fig_2_1D.Set_fig_1(Exp);
            p   = [p(1),p(2)-1];
            l_1 = true;
            l_2 = true;
            
            if ~Exp
                Fig = 3;
                figure(Fig(1)); set(gcf,'Units','pixels','Position',[150,100,1250,600]);
                subplot(1,2,1);
                Fig_2_1D.Plot_12_1(fig,l_1,1,msh,p(1),pde.e.t.f(:,1),ttm{1});
                subplot(1,2,2);
                Fig_2_1D.Plot_12_1(fig,l_2,2,msh,p(2),pde.e.t.f(:,2),ttm{2});
            else
                Fig = [3,4];
                figure(Fig(1)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_12_1(fig,l_1,1,msh,p(1),pde.e.t.f(:,1),ttm{1});
                Fig_Tools_1D.Export_PDF(F_1,D_1);
                figure(Fig(2)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_12_1(fig,l_2,2,msh,p(2),pde.e.t.f(:,2),ttm{2});
                Fig_Tools_1D.Export_PDF(F_2,D_2);
            end
        end
        % >> 1. -----------------------------------------------------------
        function [] = Plot_12_1(f,l,i,msh,p,et,ttm)
            %  > Auxiliary variables.
            C     = linspecer(9,'qualitative');
            [L{1},X] = Fig_2_1D.Set_Labels_1(i);
            trsh  = 10e-16;
            
            %  > Plot.
            hold on;
            k     = 1;
            P{1}  = plot(msh.f.Xv,et,':s','Color',C(1,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(1,:),'MarkerSize',f.MS_1);
            for j = 1:size(ttm,2)
                %  > Remove elements below a given treshold.
                ij = find(ttm(:,j)<trsh);
                if ~isempty(ij)
                    ttm(ij,j) = 0;
                end
                if ~all(ttm(:,j) < trsh)
                    k    = k+1;
                    P{k} = plot(msh.f.Xv,ttm(:,j),'--o','Color',C(j+1,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(j+1,:),'MarkerSize',f.MS_1);
                    switch i
                        case 1
                            L{k} = join(["$|\tau_{\left(",num2str(p+j),"\right)}^{\phi}|$"]);
                        case 2
                            L{k} = join(["$|\tau_{\left(",num2str(p+j),"\right)}^{\nabla\phi}|$"]);
                    end
                end
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',f.FT_1,'NumColumns',1);
            %  > Axis.
            if l
                set(gca,'YScale','log');
            end
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",X,f.FT_2,f.FT_3);
        end
        % >> 2. -----------------------------------------------------------
        %  > 2.1. ---------------------------------------------------------
        function [X,Y] = Set_Labels_1(i)
            switch i
                case 1
                    X = '$|\tau_{\phantom{(j)}}^{\phi}|$';
                    Y = "$\textrm{Error magnitude}, |\tau^{\phi}|$";
                    
                case 2
                    X = '$|\tau_{\phantom{(j)}}^{\nabla\phi}|$';
                    Y = "$\textrm{Error magnitude}, |\tau^{\nabla\phi}|$";
            end
        end
        %  > 2.2. ---------------------------------------------------------
        function [fig] = Set_fig_1(Exp)
            if ~Exp
                fig.LW_1 = 1.75;
                fig.LW_2 = 1.5;
                fig.MS_1 = 3.0;
                fig.FT_1 = 12.5;
                fig.FT_2 = 12.5;
                fig.FT_3 = 12.5;
            else
                fig.LW_1 = 3.0;
                fig.LW_2 = 2.5;
                fig.MS_1 = 5.0;
                fig.FT_1 = 25.0;
                fig.FT_2 = 30.0;
                fig.FT_3 = 30.0;
            end
        end
        %% > #2.
        function [] = WrapUp_Fig_2_2_1D(msh,pde,p,dfn_1,dfn_2,ttm_1,ttm_2)
            %  > Plot/Export.
            Exp = false;
            F_1 = "FN_1";
            F_2 = "FN_2";
            F_1 = "FN_3";
            F_2 = "FN_4";
            D_1 = "../[Figures]/[1D]/Fig_2";
            D_2 = "../[Figures]/[1D]/Fig_2";
            D_3 = "../[Figures]/[1D]/Fig_2";
            D_4 = "../[Figures]/[1D]/Fig_2";
            %  > Properties.
            fig = Fig_2_1D.Set_fig_2(Exp);
            p   = [p(1),p(2)-1];
            l_1 = true;
            l_2 = true;
            l_3 = true;
            l_4 = true;
            
            if ~Exp
                Fig = [5,6];
                figure(Fig(1)); set(gcf,'Units','pixels','Position',[150,100,1250,600]);
                subplot(1,2,1);
                Fig_2_1D.Plot_12_2(fig,l_1,1,msh,p(1),pde.e.t.f(:,1),ttm_1{1},ttm_2{1});
                subplot(1,2,2);
                Fig_2_1D.Plot_12_2(fig,l_2,2,msh,p(2),pde.e.t.f(:,2),ttm_1{2},ttm_2{2});
                figure(Fig(2)); set(gcf,'Units','pixels','Position',[150,100,1250,600]);
                subplot(1,2,1);
                Fig_2_1D.Plot_12_3(fig,l_3,msh,p(1),dfn_1{1},dfn_2{1});
                subplot(1,2,2);
                Fig_2_1D.Plot_12_3(fig,l_4,msh,p(2),dfn_1{2},dfn_2{2});
            else
                Fig = [5,6,7,8];
                figure(Fig(1)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_12_2(fig,l_1,msh,p(1),pde.e.t.f(:,1),ttm_1{1},ttm_2{1});
                Fig_Tools_1D.Export_PDF(F_1,D_1);
                figure(Fig(2)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_12_2(fig,l_2,msh,p(2),pde.e.t.f(:,2),ttm_1{2},ttm_2{2});
                Fig_Tools_1D.Export_PDF(F_2,D_2);
                figure(Fig(3)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_12_3(fig,l_3,msh,p(1),dfn_1{1},dfn_2{1});
                Fig_Tools_1D.Export_PDF(F_3,D_3);
                figure(Fig(4)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_12_3(fig,l_4,msh,p(2),dfn_1{2},dfn_2{2});
                Fig_Tools_1D.Export_PDF(F_4,D_4);
            end
        end
        % >> 1. -----------------------------------------------------------
        %  > 1.1. ---------------------------------------------------------
        function [] = Plot_12_2(f,l,i,msh,p,et,ttm_1,ttm_2)
            %  > Auxiliary variables.
            C     = linspecer(9,'qualitative');
            [X,Y] = Fig_2_1D.Set_Labels_1(i);
            trsh  = 10e-12;
            nsh   = 10;
            
            %  > Plot.
            hold on;
            k   = 0;
            im  = abs(et) > trsh;
            PX  = plot(msh.f.Xv(im),et(im),':s','Color',C(1,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(1,:),'MarkerSize',f.MS_1);
            for j = 1:size(ttm_1,2)
                %  > Remove elements below a given treshold.
                ij = abs(ttm_1(:,j)) > trsh;
                if length(find(ttm_1(:,j) > trsh)) > nsh
                    k    = k+1;
                    P{k} = plot(msh.f.Xv(ij),ttm_1(ij,j),'-o','Color',C(j+1,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(j+1,:),'MarkerSize',f.MS_1);
                    switch i
                        case 1
                            L{k} = join(["$|\tau_{\left(",num2str(p+j),"\right)}^{\phi\phantom{^{(p)}}}|$"]);
                        case 2
                            L{k} = join(["$|\tau_{\left(",num2str(p+j),"\right)}^{\nabla\phi}|$"]);
                    end
                end
            end
            for j = 1:size(ttm_2,2)
                %  > Remove elements below a given treshold.
                ik = abs(ttm_2(:,j)) > trsh;
                if length(find(ttm_2(:,j) > trsh)) > nsh
                    k    = k+1;
                    P{k} = plot(msh.f.Xv(ik),ttm_2(ik,j),'--v','Color',C(j+1,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(j+1,:),'MarkerSize',f.MS_1);
                    switch i
                        case 1
                            L{k} = join(["$|\tau_{\left(",num2str(p+j),"\right)}^{\phi^{(p)}}|$"]);
                        case 2
                            L{k} = join(["$|\tau_{\left(",num2str(p+j),"\right)}^{\nabla\phi^{(p)}}|$"]);
                    end
                end
            end
            LX{1,1} = X;
            switch i
                case 1
                    bbox = [-0.0070,-0.0600];
                case 2
                    bbox = [-0.0070,-0.0975];
            end
            legendflex(PX,LX,'buffer',bbox,'bufferunit', 'normalized',...
                'Interpreter','latex','Location','Northeast','FontSize',f.FT_1,'NumColumns',1);          
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',f.FT_1,'NumColumns',2);
            %  > Axis.
            if l
                set(gca,'YScale','log');
            end
            ylim([trsh,10^0]);
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",Y,f.FT_2,f.FT_3);
        end
        %  > 1.2. ---------------------------------------------------------
        function [] = Plot_12_3(f,l,msh,p,dfn_1,dfn_2)
            %  > Auxiliary variables.
            C    = linspecer(9,'qualitative');
            X    = "$\textrm{Derivative magnitude}, |(\underbrace{\nabla\nabla\ldots\nabla}_{p-1}\phi)_{f}|$";
            trsh = 10e-16;
            
            %  > Plot.
            hold on;
            k     = 0;
            for j = 1:size(dfn_1,2)
                %  > Remove elements below a given treshold.
                ij = abs(dfn_1(:,j)) > trsh;
                if ~all(dfn_1(:,j) < trsh)
                    k    = k+1;
                    P{k} = plot(msh.f.Xv(ij),abs(dfn_1(ij,j)),'-s','Color',C(k,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(j+1,:),'MarkerSize',f.MS_1);
                    L{k} = join(["$|\nabla^{\left(",num2str(p+j),"\right)}\phi^{\phantom{\left(p\right)}}|$"]);
                end
            end
            for j = 1:size(dfn_2,2)
                %  > Remove elements below a given treshold.
                ij = abs(dfn_2(:,j)) > trsh;
                if ~all(dfn_2(:,j) < trsh)
                    k    = k+1;
                    P{k} = plot(msh.f.Xv(ij),abs(dfn_2(ij,j)),'--o','Color',C(k,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(j+1,:),'MarkerSize',f.MS_1);
                    L{k} = join(["$|\nabla^{\left(",num2str(p+j),"\right)}\phi^{\left(p\right)}|$"]);
                end
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',f.FT_1,'NumColumns',1);
            %  > Axis.
            if l
                set(gca,'YScale','log');
            end
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",X,f.FT_2,f.FT_3);
        end
        % >> 2.
        % >> 2. -----------------------------------------------------------
        %  > 2.1. ---------------------------------------------------------
        function [fig] = Set_fig_2(Exp)
            if ~Exp
                fig.LW_1 = 1.75;
                fig.LW_2 = 1.5;
                fig.MS_1 = 3.0;
                fig.FT_1 = 12.5;
                fig.FT_2 = 12.5;
                fig.FT_3 = 12.5;
            else
                fig.LW_1 = 3.0;
                fig.LW_2 = 2.5;
                fig.MS_1 = 5.0;
                fig.FT_1 = 25.0;
                fig.FT_2 = 30.0;
                fig.FT_3 = 30.0;
            end
        end
    end
end