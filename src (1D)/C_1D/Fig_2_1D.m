classdef Fig_2_1D
    methods (Static)
        %% > #1.
        function [] = WrapUp_Fig_2_1_1D(msh,et,tm,p)
            %  > Plot/Export.
            Exp = 0;
            F_1 = "FN_1";
            F_2 = "FN_2";
            D_1 = "../[Figures]/[1D]/Fig_2";
            D_2 = "../[Figures]/[1D]/Fig_2";
            %  > Properties.
            Fig = [1,2];
            fig = Fig_2_1D.Set_fig_1(Exp);
            p   = [p(1),p(2)-1];
            l_1 = 1;
            l_2 = 1;
            
            if ~Exp
                %  > 1/2.
                figure(Fig(1)); set(gcf,'Units','pixels','Position',[150,100,1250,600]);
                subplot(1,2,1);
                Fig_2_1D.Plot_1(fig,l_1,1,msh,p(1),et.f_abs(:,1),tm{1});
                subplot(1,2,2);
                Fig_2_1D.Plot_1(fig,l_2,2,msh,p(2),et.f_abs(:,2),tm{2});
            else
                %  > 1.
                figure(Fig(1)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_1(fig,l_1,1,msh,p(1),et.f_abs(:,1),tm{1});
                Fig_Tools_1D.Export_PDF(F_1,D_1);
                %  > 2.
                figure(Fig(2)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_1(fig,l_2,2,msh,p(2),et.f_abs(:,2),tm{2});
                Fig_Tools_1D.Export_PDF(F_2,D_2);
            end
        end
        % >> 1. -----------------------------------------------------------
        function [] = Plot_1(f,l,i,msh,p,et,tm)
            %  > Auxiliary variables.
            C        = linspecer(9,'qualitative');
            [L{1},X] = Fig_2_1D.Set_Labels_1(i);
            trsh     = 10e-12;
            
            hold on;
            k     = 1;
            P{1}  = plot(msh.f.Xv,et,':s','Color',C(1,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(1,:),'MarkerSize',f.MS_1);
            for j = 1:size(tm,2)
                %  > tm.
                m = find(tm(:,j) > trsh);
                if ~all (tm(:,j) < trsh)
                    k    = k+1;
                    P{k} = plot(msh.f.Xv,tm(:,j),':o','Color',C(j+1,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(j+1,:),'MarkerSize',f.MS_1);
                    switch i
                        case 1
                            L{k} = join(["$|\tau_{f}^{\phi^{\left(",num2str(p+j),"\right)}}|$"]);
                        case 2
                            L{k} = join(["$|\tau_{f}^{\nabla\phi^{\left(",num2str(p+j),"\right)}}|$"]);
                    end
                    if l
                        y_min(k-1) = min(tm(m,j));
                        y_max(k-1) = max(tm(m,j));
                    end
                end
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',f.FT_1,'NumColumns',1);
            %  > Axis.
            if l
                set(gca,'YScale','log');
                ylim([10.^(ceil(log10(min(y_min)))-1),...
                      10.^(ceil(log10(max(y_max)))+1)]);
            end
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",X,f.FT_2,f.FT_3);
        end
        % >> 2. -----------------------------------------------------------
        %  > 2.1. ---------------------------------------------------------
        function [X,Y] = Set_Labels_1(i)
            switch i
                case 1
                    X = '$|\tau_{f}^{\phi^{\phantom{\left(j\right)}}}|$';
                    Y = "$\textrm{Error magnitude}, |\tau_{f}^{\phi}|$";
                    
                case 2
                    X = '$|\tau_{f}^{\nabla\phi^{\phantom{\left(j\right)}}}|$';
                    Y = "$\textrm{Error magnitude}, |\tau_{f}^{\nabla\phi}|$";
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
                fig.FT_1 = 21.5;
                fig.FT_2 = 25.0;
                fig.FT_3 = 25.0;
            end
        end
        %% > #2.
        function [] = WrapUp_Fig_2_2_1D(msh,p,df,tm)
            %  > Plot/Export.
            Exp = false;
            F_1 = "FN_1";
            F_2 = "FN_2";
            F_3 = "FN_3";
            F_4 = "FN_4";
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
                Fig = [1,2];
                figure(Fig(1)); set(gcf,'Units','pixels','Position',[150,100,1250,600]);
                subplot(1,2,1);
                Fig_2_1D.Plot_2_1(fig,l_1,1,msh,p(1),tm.a{1},tm.n{1});
                subplot(1,2,2);
                Fig_2_1D.Plot_2_1(fig,l_2,2,msh,p(2),tm.a{2},tm.n{2});
                figure(Fig(2)); set(gcf,'Units','pixels','Position',[150,100,1250,600]);
                subplot(1,2,1);
                Fig_2_1D.Plot_2_2(fig,l_3,msh,p(1),df.a{1},df.n{1});
                subplot(1,2,2);
                Fig_2_1D.Plot_2_2(fig,l_4,msh,p(2),df.a{2},df.n{2});
            else
                Fig = [1,2,3,4];
                figure(Fig(1)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_2_1(fig,l_1,msh,p(1),tm.a{1},tm.n{1});
                Fig_Tools_1D.Export_PDF(F_1,D_1);
                figure(Fig(2)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_2_1(fig,l_2,msh,p(2),tm.a{2},tm.n{2});
                Fig_Tools_1D.Export_PDF(F_2,D_2);
                figure(Fig(3)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_2_2(fig,l_3,msh,p(1),df.a{1},df.n{1});
                Fig_Tools_1D.Export_PDF(F_3,D_3);
                figure(Fig(4)); set(gcf,'Units','pixels','Position',[250,100,1050,650]);
                Fig_2_1D.Plot_2_2(fig,l_4,msh,p(2),df.a{2},df.n{2});
                Fig_Tools_1D.Export_PDF(F_4,D_4);
            end
        end
        % >> 1. -----------------------------------------------------------
        %  > 1.1. ---------------------------------------------------------
        function [] = Plot_2_1(f,l,i,msh,p,tm_a,tm_n)
            %  > Auxiliary variables.
            C     = linspecer(9,'qualitative');
            [~,Y] = Fig_2_1D.Set_Labels_1(i);
            trsh  = 10e-12;
            nsh   = 10;
            
            %  > Plot.
            hold on;
            k  = 0;
            for j = 1:size(tm_a,2)
                %  > tm_a.
                m = abs(tm_a(:,j)) > trsh;
                if length(find(tm_a(:,j) > trsh)) > nsh
                    k    = k+1;
                    P{k} = plot(msh.f.Xv,tm_a(:,j),'-o','Color',C(j+1,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(j+1,:),'MarkerSize',f.MS_1);
                    switch i
                        case 1
                            L{k} = join(["$|\tau_{\left(",num2str(p+j),"\right)}^{\phi\phantom{^{(p)}}}|$"]);
                        case 2
                            L{k} = join(["$|\tau_{\left(",num2str(p+j),"\right)}^{\nabla\phi}|$"]);
                    end
                    if l
                        y_min(k) = min(tm_a(m,j));
                        y_max(k) = max(tm_a(m,j));
                    end
                end
            end
            for j = 1:size(tm_n,2)
                %  > tm_n.
                n = abs(tm_n(:,j)) > trsh;
                if length(find(tm_n(:,j) > trsh)) > nsh
                    k    = k+1;
                    P{k} = plot(msh.f.Xv,tm_n(:,j),'--v','Color',C(j+1,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(j+1,:),'MarkerSize',f.MS_1);
                    switch i
                        case 1
                            L{k} = join(["$|\tau_{\left(",num2str(p+j),"\right)}^{\phi^{(p)}}|$"]);
                        case 2
                            L{k} = join(["$|\tau_{\left(",num2str(p+j),"\right)}^{\nabla\phi^{(p)}}|$"]);
                    end
                    if l
                        y_min(k) = min(tm_n(n,j));
                        y_max(k) = max(tm_n(n,j));
                    end
                end
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',f.FT_1,'NumColumns',2);
            %  > Axis.
            if l
                set(gca,'YScale','log');
                ylim([10.^(ceil(log10(min(y_min)))-1),...
                      10.^(ceil(log10(max(y_max)))+2)]);
            end
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",Y,f.FT_2,f.FT_3);
        end
        %  > 1.2. ---------------------------------------------------------
        function [] = Plot_2_2(f,l,msh,p,dfa,dfn)
            %  > Auxiliary variables.
            C    = linspecer(9,'qualitative');
            X    = "$\textrm{Derivative magnitude}$";
            trsh = 10e-10;
            
            %  > Plot.
            hold on;
            k     = 0;
            for j = 1:size(dfa,2)
                %  > dfa.
                m = abs(dfa(:,j)) > trsh;
                if ~all(dfa(:,j)  < trsh)
                    k    = k+1;
                    P{k} = plot(msh.f.Xv,abs(dfa(:,j)),'-s','Color',C(k,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(j+1,:),'MarkerSize',f.MS_1);
                    L{k} = join(["$|\nabla^{\left(",num2str(p+j),"\right)}\phi\phantom{,}|$"]);
                    if l
                        y_min(k) = min(abs(dfa(m,j)));
                        y_max(k) = max(abs(dfa(m,j)));
                    end
                end
            end
            for j = 1:size(dfn,2)
                %  > dfn.
                n = abs(dfn(:,j)) > trsh;
                if ~all(dfn(:,j)  < trsh)
                    k    = k+1;
                    P{k} = plot(msh.f.Xv,abs(dfn(:,j)),'--o','Color',C(k,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(j+1,:),'MarkerSize',f.MS_1);
                    L{k} = join(["$|\nabla^{\left(",num2str(p+j),"\right)}\phi^{\left(p\right)}|$"]);
                    if l
                        y_min(k) = min(abs(dfn(n,j)));
                        y_max(k) = max(abs(dfn(n,j)));
                    end
                end
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',f.FT_1,'NumColumns',2);
            %  > Axis.
            if l
                set(gca,'YScale','log');
                ylim([10.^(ceil(log10(min(y_min)))-1),...
                      10.^(ceil(log10(max(y_max)))+2)]);
            end
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",X,f.FT_2,f.FT_3);
        end
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
        %% > #3.
        function [] = WrapUp_Fig_2_3_1D(msh,pde,p)
            %  > Plot/Export.
            Exp  = false;
            F_1  = "2_3";
            F_2  = "2_3";
            D_1  = "../[Figures]/[1D]/Fig_2";
            D_2  = "../[Figures]/[1D]/Fig_2";
            %  > Properties.
            fig  = Fig_Tools_1D.Set_fig(Exp);
            plot = [0,1];

            if ~Exp
                Fig = [1,2,3];
                if plot(1)
                    figure(Fig(1)); set(gcf,'Units','pixels','Position',[150,100,1250,600]);
                    subplot(1,2,1);
                    Fig_2_1D.Plot_3_1(fig,1,msh,pde,p);
                    subplot(1,2,2);
                    Fig_2_1D.Plot_3_1(fig,2,msh,pde,p);
                end
                if plot(2)
                    figure(Fig(2)); set(gcf,'Units','pixels','Position',[150,100,1250,600]);
                    subplot(1,2,1);
                    Fig_2_1D.Plot_3_2(fig,1,msh,pde,p);
                    subplot(1,2,2);
                    Fig_2_1D.Plot_3_2(fig,2,msh,pde,p);
                end
            else
                Fig = [1,2,3];
            end
        end
        % >> 1. -----------------------------------------------------------
        %  > 1.1. ---------------------------------------------------------
        function [] = Plot_3_1(f,i,msh,pde,p)
            %  > Auxiliary variables.
            C    = linspecer(9,'qualitative');
            switch i
                case 1
                    str = "\phi";
                case 2
                    str = "\nabla\phi";
                otherwise
                    return;
            end
            X    = join(["$\textrm{Error magnitude}, |\tau_{f}^{",str,"}|$"]);
            trsh = 10e-12;
            
            %  > Tau_f(m/n): j-th polynomial approximation w/ l-th order solution.
            hold on;
            [m,n] = size(pde.et.df.x); k = 0; q = 0;
            for j = 1:m
                q = q+1;
                for l = 1:m
                    o           = pde.et.fv.e{j,l}.f(:,i) > trsh; pde.et.fv.e{j,l}.f(~o,i) = 0;
                    k           = k+1; 
                    M           = Fig_Tools_1D.M(q);
                    P1     {k}  = plot(msh.f.Xv,pde.et.fv.e{j,l}.f(:,i),M,'Color',C(q,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(q,:),'MarkerSize',f.MS_1);
                    L1   {1,k}  = Fig_2_1D.Set_Labels_2_3_1(j,l,p,str);
                    L1   {1,k} = convertStringsToChars(L1{1,k});
                    [A(k),B(k)] = MinMaxElem(pde.et.fv.e{j,l}.f(o,i),'finite');
                end
            end
            %  > Tau_f(m/n-m/m).
            k = 0;
            for j = 1:m
                q = q+1;
                N = setdiff(p,p(j));
                for l = 1:n
                    o           = pde.et.df.x{j,l}(:,i) > trsh; pde.et.df.x{j,l}(~o,i) = 0;
                    k           = k+1;
                    M           = Fig_Tools_1D.M(q);
                    P2      {k} = plot(msh.f.Xv,pde.et.df.x{j,l}(:,i),M,'Color',C(q,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(q,:),'MarkerSize',f.MS_1);
                    L2    {1,k} = Fig_2_1D.Set_Labels_2_3_2(p(j),p(j),p(j),N(l),str);
                    L2    {1,k} = convertStringsToChars(L2{1,k});
                    [A(k),B(k)] = MinMaxElem(pde.et.df.x{j,l}(o,i),'finite');
                end 
            end
            switch i
                case 1
                    bbox_1 = [-0.090,-0.025];
                    bbox_2 = [-0.010,-0.025];
                case 2
                    bbox_1 = [-0.097,-0.025];
                    bbox_2 = [-0.010,-0.025];
                otherwise
                    return;
            end
            legendflex([P1{:}],L1,'buffer',bbox_1,'bufferunit','normalized',...
                'Interpreter','latex','Location','Northeast','FontSize',f.FT_1,'NumColumns',1);
            legendflex([P2{:}],L2,'buffer',bbox_2,'bufferunit', 'normalized',...
                'Interpreter','latex','Location','Northeast','FontSize',f.FT_1,'NumColumns',1);
            %  > Axis.
            set(gca,'YScale','log');
            ylim([10.^(ceil(log10(min(A)))-1),...
                  10.^(ceil(log10(max(B)))+1)]);
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",X,f.FT_2,f.FT_3);
        end
        %  > 1.2. ---------------------------------------------------------
        function [] = Plot_3_2(f,i,msh,pde,p)
            %  > Auxiliary variables.
            C    = linspecer(9,'qualitative');
            switch i
                case 1
                    str = "\phi";
                case 2
                    str = "\nabla\phi";
                otherwise
                    return;
            end
            X    = join(["$\textrm{Error magnitude}, |\tau_{f}^{",str,"}|$"]);
            trsh = 10e-12;
            
            %  > Tau_f(m).
            hold on;
            [m,n] = size(pde.et.df.x); 
            
            %  > Tau_f(m/n-m/m).
            k = 0; q = 1;
            for j = 1:n
                o           = pde.et.av{j}.f_abs(:,i) > trsh; pde.et.av{j}.f_abs(~o,i) = 0;
                k           = k+1;
                M           = Fig_Tools_1D.M(q);
                P{k}        = plot(msh.f.Xv,pde.et.av{j}.f_abs(:,i),M,'Color',C(q,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(q,:),'MarkerSize',f.MS_1);
                L{k}        = Fig_2_1D.Set_Labels_2_3_3(p(j),str);
                [A(k),B(k)] = MinMaxElem(pde.et.av{j}.f_abs(o,i),'finite');
            end
            %  > Tau_f(m-n).
            q = q+1;
            for j = 1:n
                o           = pde.et.df.a{j}(:,i) > trsh; pde.et.df.a{j}(~o,i) = 0;
                k           = k+1;
                M           = Fig_Tools_1D.M(q);
                P{k}        = plot(msh.f.Xv,pde.et.df.a{j}(:,i),M,'Color',C(q,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(5,:),'MarkerSize',f.MS_1);
                L{k}        = Fig_2_1D.Set_Labels_2_3_4(p(j),p(j+1),str);
                [A(k),B(k)] = MinMaxElem(pde.et.df.a{j}(o,i),'finite');
            end 
            for j = 1:m
                q = q+1;
                N = setdiff(p,p(j));
                for l = 1:n
                    o           = pde.et.df.x{j,l}(:,i) > trsh; pde.et.df.x{j,l}(~o,i) = 0;
                    k           = k+1;
                    M           = Fig_Tools_1D.M(q);
                    P      {k}  = plot(msh.f.Xv,pde.et.df.x{j,l}(:,i),M,'Color',C(q,:),'LineWidth',f.LW_1,'MarkerFaceColor',C(q,:),'MarkerSize',f.MS_1);
                    L      {k}  = Fig_2_1D.Set_Labels_2_3_2(p(j),p(j),p(j),N(l),str);
                    [A(k),B(k)] = MinMaxElem(pde.et.df.x{j,l}(o,i),'finite');
                end 
            end
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',f.FT_1,'NumColumns',1);
            %  > Axis.
            set(gca,'YScale','log');
            ylim([10.^(ceil(log10(min(A)))-1),...
                10.^(ceil(log10(max(B)))+1)]);
            Fig_Tools_1D.ChangeLook_1D(true,true,true,msh.f.Xv,10,"$x$",X,f.FT_2,f.FT_3);
        end
        % >> 2. -----------------------------------------------------------
        %  > 2.1. ---------------------------------------------------------
        %  > 2.1.1. -------------------------------------------------------
        function [L] = Set_Labels_2_3_1(a,b,p,str)
            L = join(["$|\tau_{f}^{",str,"^{\left(",num2str(p(a)),"/",num2str(p(b)),"\right)}}|$"]);
        end
        %  > 2.1.2. -------------------------------------------------------
        function [L] = Set_Labels_2_3_2(a,b,c,d,str)
            L = join(["$|\tau_{f}^{",str,"^{\left(",num2str(a),"/",num2str(b),"-",num2str(c),"/",num2str(d),"\right)}}|$"]);
        end
        %  > 2.1.3. -------------------------------------------------------
        function [L] = Set_Labels_2_3_3(a,str)
            L = join(["$|\tau_{f}^{",str,"^{\left(",num2str(a),"\right)\phantom{-:}}}|$"]);
        end
        %  > 2.1.4. -------------------------------------------------------
        function [L] = Set_Labels_2_3_4(a,b,str)
            L = join(["$|\tau_{f}^{",str,"^{\left(",num2str(a),"-",num2str(b),"\right)}}|$"]);
        end
    end
end