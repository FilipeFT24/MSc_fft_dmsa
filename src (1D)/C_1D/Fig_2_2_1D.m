classdef Fig_2_2_1D
    methods (Static)
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