classdef Fig_Tools
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [fig] = Set_fig(opt,x)
            fig.x            = x;
            fig.C            = linspecer(9,'qualitative');       %  > Colors.
            fig.LW           =  2.25;                            %  > Line width.
            fig.MS           =  8.50;                            %  > Marker size.
            fig.FA           =  0.25;                            %  > FaceAlpha.
            if ~opt
                fig.FS{1}(1) = 35.00;                            %  > x-label.        (*)
                fig.FS{1}(2) = 35.00;                            %  > y-label.        (*)
            else
                fig.FS{1}(1) = 22.50;                            %  > x-label.        (*)
                fig.FS{1}(2) = 22.50;                            %  > y-label.        (*)
            end
            fig.FS    {2}    = 15.00;                            %  > Axis.           (*)
            fig.FS    {3}    = 20.00;                            %  > Legend.
            fig.FS    {4}    = 13.50;                            %  > Colormap label. (*)
            fig.FS    {5}    = 20.00;                            %  > Colormap title. (*)
            fig.TrshV        = 1E-16;                            %  > Error treshold value.
            switch x.l
                case 1, fig.L{1} = "$x$";                        %  > x-axis label.
                        fig.L{2} = "$\textrm{Error magnitude}$"; %  > y-axis label.
                case 2, fig.L{1} = "$\textrm{NNZ}$";             %  > x-axis label.
                        fig.L{2} = "$\textrm{Error magnitude}$"; %  > y-axis label.
                case 3, fig.L{1} = "$x$";                        %  > x-axis label.
                        fig.L{2} = "$y$";                        %  > y-axis label.
                otherwise
                    return;
            end
            %  > Markers/line styles.
            %  M = ['+','o','*','x','v','d','^','s','>','<'];
            %  L = ['-','--',':','-.'];
        end
        % >> 1.2. ---------------------------------------------------------
        %  > 1.2.1. -------------------------------------------------------
        function [S] = S1(i)
            switch i
                case 1, S = "\phi\scriptscriptstyle{,C}";
                case 2, S = "\phi\scriptscriptstyle{,D}";
                case 3, S = "\phi";
                otherwise
                    return;
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        function [S] = S2(i)
            switch i
                case 1, S = "1";
                case 2, S = "2";
                case 3, S = "\infty";
                otherwise
                    return;
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 2.1.1. -------------------------------------------------------
        function [] = Map_1D(fig,X,Z)
            %  > Axes.
            [XY.Lim(1,1),XY.Lim(1,2)] = MinMaxElem(X);
            [XY.Lim(2,1),XY.Lim(2,2)] = ...
                deal(10.^(ceil(log10(Z.lim(1)))-Z.DY(1)),10.^(ceil(log10(Z.lim(2)))+Z.DY(2)));
            set(get(gca,'XAxis'),'Fontsize',fig.FS{2}); set(gca,'box','on'); xlabel(fig.L{1},'Fontsize',fig.FS{1}(1));
            set(get(gca,'YAxis'),'Fontsize',fig.FS{2}); set(gca,'box','on'); ylabel(fig.L{2},'Fontsize',fig.FS{1}(2));
            if ~Z.Plot
                set(gca,'box','on','XLim',XY.Lim(1,:),'XTick',linspace(XY.Lim(1,1),XY.Lim(1,2),11));
            else
                if numel(unique(XY.Lim(1,:))) ~= 1
                    set(get(gca,'XAxis'),'Exponent',3); set(gca,'XLim',XY.Lim(1,:),'XScale','Log');
                end
                %  > Grid.
                grid(gca,'on');
                set (gca,'GridLineStyle','-','MinorGridLinestyle',':');
            end
            set(gca,'YLim',XY.Lim(2,:),'YScale','Log');
            %  > Figure.
            set(gcf,'Color','w','Windowstate','Maximized');
            %  > Interpreter.
            set(0,'DefaultAxesTickLabelInterpreter','Latex');
            set(0,'DefaultLegendInterpreter','Latex');
            %  > Legend.
            legend([Z.P{:}],[Z.L{:}],'Location','Northeast','FontSize',fig.FS{3},'NumColumns',Z.NC);
            %  > Maximize...
            pbaspect([1,1,1]);
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        function [LY,P,Y_lim] = Var_1D_1(fig,m,LX,X,Y)
            %  > Auxiliary variables.
            a = size(X,2);
            b = size(Y,2);
            if a < b
                X = repelem(X,1,b);
            end
            
            k = 0;
            hold on;
            for i = 1:b
                j{i} = Y(:,i) > 0; Y(~j{i},i) = 0;
                if nnz(j{i})  > fig.TrshV
                    k                 = k+1;
                    LY{k}             = LX{i};
                    P {k}             = plot(X(:,i),Y(:,i),m(i),'Color',fig.C(k,:),'LineWidth',fig.LW,'MarkerFaceColor',fig.C(k,:),'MarkerSize',fig.MS);
                    [YV(k,1),YV(k,2)] = MinMaxElem(Y(j{i},i));
                end
            end
            [Y_lim(1),Y_lim(2)] = MinMaxElem(YV);
        end
        %  > 2.2.2. -------------------------------------------------------
        function [LY,P1,YV] = Var_1D_2(fig,M1,LX,X1,Y1)
            %  > Auxiliary variables.
            [X2(1),X2(2)] = MinMaxElem(X1);
            
            k = 0;
            hold on;
            for i = 1:size(Y1,2)
                k                 = k+1;
                LY{k}             = LX{i};
                P1{k}             = fplot(Y1{1,i},X2,M1(1,i),'Color',fig.C(k,:),'LineWidth',fig.LW,'MarkerFaceColor',fig.C(k,:),'MarkerSize',fig.MS);
                P2{k}             =  plot(X1,Y1{2,i},M1(2,i),'Color',fig.C(k,:),'LineWidth',fig.LW,'MarkerFaceColor',fig.C(k,:),'MarkerSize',fig.MS);
                [YV(k,1),YV(k,2)] = MinMaxElem(Y1{2,i});
            end
        end
        %  > 2.2.3. -------------------------------------------------------
        function [LY,P,YV] = Var_1D_3(fig,M,LX,X,Y)
            k = 0;
            hold on;
            for i = 1:size(Y,2)
                j = Y(i) > fig.trsh; Y(~j) = 0;
                if j
                    k                 = k+1;
                    P {k}             = line([X(1),X(end)],[Y(i),Y(i)],'Color',fig.C(i,:),'Linewidth',fig.LW,'Linestyle',M(i));
                    LY{k}             = LX{i};
                    [YV(k,1),YV(k,2)] = MinMaxElem(Y(i));
                end
            end
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        function [] = Map_2D(fig,msh,Z)
            %  > Axes.
            BBOX  = get(gca,'Position');
            XY.V  = cat(1,msh.f.xy.v{:});
            for i = 1:size(XY.V,2)
                [XY.Lim(i,1),XY.Lim(i,2)] = MinMaxElem(XY.V(:,i));
            end
            set(get(gca,'XAxis'),'Fontsize',fig.FS{2}); set(gca,'box','on','XLim',XY.Lim(1,:),'XTick',[]); xlabel(fig.L{1},'Fontsize',fig.FS{1}(1));
            set(get(gca,'YAxis'),'Fontsize',fig.FS{2}); set(gca,'box','on','YLim',XY.Lim(2,:),'YTick',[]); ylabel(fig.L{2},'Fontsize',fig.FS{1}(2));
            %  > Colormap.
            if exist('P','Var')
                colorbar (gca,'box','on','Fontsize',fig.FS{4},'Position',[BBOX(1)+BBOX(3)+0.0075,0.1755,0.0225,0.6855]); AdvancedColormap('thermal');
                set      (gca,'ColorScale','Log');
                caxis    (Z.Lim);
                if ~Z.e
                    title(Z.Title,'Fontsize',fig.FS{5});
                end
            end
            %  > Figure.
            set(gcf,'Color','w','Windowstate','Maximized');
            %  > Grid.
            Fig_Tools.Var_2D_1(msh.f.xy.v,"k",1E-3,"-");
            %  > Interpreter.
            set(0,'DefaultAxesTickLabelInterpreter','Latex');
            set(0,'DefaultLegendInterpreter','Latex');
            %  > Maximize...
            pbaspect([1,1,1]);
        end
        % >> 3.2. ---------------------------------------------------------
        %  > 3.2.1. -------------------------------------------------------
        function [] = Var_2D_1(V,C,LW,M)
            hold on;
            for i = 1:size(V,1)
                plot(V{i}(:,1),V{i}(:,2),M,'Color',C,'LineWidth',LW);
            end
        end
        %  > 3.2.2. -------------------------------------------------------
        function [] = Var_2D_2(V,C,M,MS)
            if numel(MS) ~= size(V,1)
                MS = ones(size(V,1),1).*MS;
            end
            hold on;
            for i = 1:size(V,1)
                plot(V(:,1),V(:,2),M,'Color',C,'MarkerFaceColor',C,'MarkerSize',MS(i));
            end
        end
        %  > 3.2.3. -------------------------------------------------------
        function [] = Var_2D_3(V,C,LW,M,MS)
            if numel(MS) ~= size(V,1)
                MS = ones(size(V,1),1).*MS;
            end
            hold on;
            for i = 1:size(V,1)
                plot(V(:,1),V(:,2),M,'Color',C,'LineWidth',LW,'MarkerSize',MS(i));
            end
        end
        %  > 3.2.4. -------------------------------------------------------
        function [] = Var_2D_4(V1,V2,FA)
            hold on;
            [sz(1),sz(2)] = deal(size(V1,1),size(V2,1));
            if sz(1) > sz(2)
                V2 = repelem(V2,sz(1),1);
            end
            for i = 1:sz(1)
                patch(V1{i}(:,1),V1{i}(:,2),V2(i,:),'FaceAlpha',FA,'Linestyle','None');
            end
        end
        %  > 3.2.5. -------------------------------------------------------
        function [] = Var_2D_5(V1,V2,C,Length)
            hold on;
            for i = 1:size(V1,1)
                for j = 1:size(V1{i},1)
                    quiver(V1{i}(j,1),V1{i}(j,2),V2{i}(j,1)./Length,V2{i}(j,2)./Length,'AutoScale','off','Color',C);
                end
            end
        end
        %  > 3.2.6. -------------------------------------------------------
        function [] = Var_2D_6(V1,V2,V3,V4,sc,sf,x,w,C,FA,LW,M,MS)
            %  > "c".
            Fig_Tools.Var_2D_2(V1,"k",M(1),MS(1));
            for i = 1:numel(sc)
                Fig_Tools.Var_2D_2(V1(sc{i},:),C(i,:),M(2),MS(2));
                Fig_Tools.Var_2D_4(V2(sc{i})  ,C(i,:)     ,FA);
            end
            %  > "f".
            Fig_Tools.Var_2D_1(V3,C(1,:),LW,"-");
            Fig_Tools.Var_2D_2(x ,C(1,:),M(1),MS(2).*w);
            for i = 1:numel(sf)
                Fig_Tools.Var_2D_2(V4(sf{i},:),C(i,:),M(2),MS(2));
            end
        end
        %  > 3.2.7. -------------------------------------------------------
        function [MS] = Convert_MS(Length)
            A         = gca;
            B         = A.Units;
            A.Units   = 'Points';
            C         = A.Position(3:4);
            A.Units   = B;
            MS_V      = C.*Length./[range(xlim(A)),range(ylim(A))];
            MS        = MS_V(2);
        end
    end
end