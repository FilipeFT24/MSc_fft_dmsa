classdef Fig_Tools
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [fig] = Set_fig(opt,x)
            fig.x            = x;
            fig.C            = linspecer(9,'qualitative');       %  > Colors.
            fig.lw           =  1.75;                            %  > Line width.
            fig.ms           =  5.00;                            %  > Marker size.
            fig.fa           =  0.25;                            %  > FaceAlpha.
            if ~opt
                fig.fs{1}(1) = 30.00;                            %  > x-label.
                fig.fs{1}(2) = 30.00;                            %  > y-label.
            else
                fig.fs{1}(1) = 22.50;                            %  > x-label.
                fig.fs{1}(2) = 22.50;                            %  > y-label.
            end
            fig.fs    {2}    = 15.00;                            %  > Axis.
            fig.fs    {3}    = 15.00;                            %  > Legend.
            fig.fs    {4}    = 11.50;                            %  > Colormap label.
            fig.fs    {5}    = 20.00;                            %  > Colormap title.
            fig.pos          = [150,100,1250,600];               %  > Position.
            fig.dir          = '[Post-processing]/[.pdf Files]'; %  > Destination folder.
            fig.fn           = 'Times New Roman';                %  > Font name.
            fig.trsh.n       = 0;                                %  > Error treshold: number of elements.
            fig.trsh.v       = 1.0e-16;                          %  > Error treshold: value.
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
        function [str] = Set_str_1(i)
            switch i
                case 1
                    str = "\phi\scriptscriptstyle{,C}";
                case 2
                    str = "\phi\scriptscriptstyle{,D}";
                case 3
                    str = "\phi";
                otherwise
                    return;
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        function [str] = Set_str_2(i)
            switch i
                case 1
                    str = "1";
                case 2
                    str = "2";
                case 3
                    str = "\infty";
                otherwise
                    return;
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 2.1.1. -------------------------------------------------------
        function [] = Map_1D(fig,x,y)
            %  > Other parameters.
            a = gca;
            f = gcf;
            box on;
            set(a,'Fontname',fig.fn)
            set(a,'Fontsize',fig.fs{2});
            set(a,'TickLabelInterpreter','Latex');
            set(a,'YScale','log');
            set(f,'Color','w');
            set(f,'Windowstate','maximized');
            %  > Axis.
            [X_lim(1),X_lim(2)] = MinMaxElem(x);
            if ~y.plot
                n  = 10;
                dX = diff(X_lim)./n; X_label = X_lim(1)-dX:dX:X_lim(2)+dX;
                set(a,'XLim',X_lim,'XTick',X_label);
            else
                grid on; 
                grid minor;
                if numel(unique(X_lim)) ~= 1
                    set(a,'XLim',X_lim);
                end
                set(a,'XScale','log'); a.XAxis.Exponent = 3;
                set(a,'GridLineStyle','-','minorgridlinestyle',':');
            end
            Y_lim (1) = 10.^(ceil(log10(y.lim(1)))-y.dY(1));
            Y_lim (2) = 10.^(ceil(log10(y.lim(2)))+y.dY(2)); ylim(Y_lim);
            xlabel(fig.L{1},'Fontsize',fig.fs{1}(1),'Interpreter','latex');
            ylabel(fig.L{2},'Fontsize',fig.fs{1}(2),'Interpreter','latex');
            %  > Legend.
            legend([y.P{:}],[y.L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',fig.fs{3},'NumColumns',y.nc);
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
                j{i} = Y(:,i) > fig.trsh.n; Y(~j{i},i) = 0;
                if nnz(j{i})  > fig.trsh.v
                    k                 = k+1;
                    LY{k}             = LX{i};
                    P {k}             = plot(X(:,i),Y(:,i),m(i),'Color',fig.C(k,:),'LineWidth',fig.lw,'MarkerFaceColor',fig.C(k,:),'MarkerSize',fig.ms);
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
        function [] = Map_2D(fig,msh,y)
            %  > Other parameters.
            a = gca;
            f = gcf;
            axis equal;
            set(a,'Fontname',fig.fn)
            set(a,'Fontsize',fig.fs{2});
            set(a,'Layer','Top')
            set(a,'TickLabelinterpreter','Latex');
            set(f,'Color','w');
            set(f,'Windowstate','maximized');
            pbaspect([1,1,1]);
            %  > Axis.
            n = [10,10];
            X = cat(1,msh.f.xy.v{:});
            for i = 1:size(X,2)
                [X_lim(i,1),X_lim(i,2)] = MinMaxElem(X(:,i)); dX(i) = diff(X_lim(i,:))./n(i); X_label(i,:) = X_lim(i,1)-dX(i):dX(i):X_lim(i,2)+dX(i);
            end
            set(a,'XLim',X_lim(1,:),'XTick',X_label(1,:));
            a.XRuler.Axle.Visible = 'off';
            a.XAxis.TickLength    = [0,0];
            set(a,'YLim',X_lim(2,:),'YTick',X_label(2,:));
            a.YRuler.Axle.Visible = 'off';
            a.YAxis.TickLength    = [0,0];
            xlabel(fig.L{1},'Fontsize',fig.fs{1}(1),'Interpreter','latex');
            ylabel(fig.L{2},'Fontsize',fig.fs{1}(2),'Interpreter','latex');
            %  > Colormap.
            if y.plot
                p                      = a.Position;
                of_p                   = 0.035;
                c                      = colorbar(a,'Location','east');
                c.FontSize             = fig.fs{4};
                c.Label.Interpreter    = 'latex';
                c.TickLabelInterpreter = 'latex';
                AdvancedColormap('thermal');
                if y.lim (1) >= fig.trsh.v
                    Y_lim(1)  = 10.^(ceil(log10(y.lim(1)))-1);
                    Y_lim(2)  = 10.^(ceil(log10(y.lim(2)))+0);
                else
                    Y_lim(1)  = fig.trsh.v;
                    Y_lim(2)  = 10.^(ceil(log10(y.lim(2)))+0);
                end
                caxis(Y_lim);
                set  (a,'colorscale','log');
                set  (c,'position',[p(1)+p(3)+of_p,0.165,of_p./2.5,0.705]);
                title(y.title,'Fontsize',fig.fs{5},'Interpreter','latex');
            end
            %  > Grid.
            Fig_Tools.Var_2D_1(msh.f.xy.v,"k","-",0.01);
        end
        % >> 3.2. ---------------------------------------------------------
        %  > 3.2.1. -------------------------------------------------------
        function [] = Var_2D_1(V,C,m,lw)
            hold on;
            for i = 1:size(V,1)
                plot(V{i}(:,1),V{i}(:,2),m,'Color',C,'LineWidth',lw);
            end
        end
        %  > 3.2.2. -------------------------------------------------------
        function [] = Var_2D_2(V,c,m,ms)
            hold on;
            if numel(ms) ~= size(V,1)
                ms = ones(size(V,1),1).*ms;
            end
            for i = 1:size(V,1)
                plot(V(:,1),V(:,2),m,'Color',c,'MarkerFaceColor',c,'MarkerSize',ms(i));
            end
        end
        %  > 3.2.3. -------------------------------------------------------
        function [] = Var_2D_3(V,c,lw,m,ms)
            hold on;
            if numel(ms) ~= size(V,1)
                ms = ones(size(V,1),1).*ms;
            end
            for i = 1:size(V,1)
                plot(V(:,1),V(:,2),m,'Color',c,'LineWidth',lw,'MarkerSize',ms(i));
            end
        end
        %  > 3.2.4. -------------------------------------------------------
        function [] = Var_2D_4(V1,V2,fa)
            hold on;
            [sz(1),sz(2)] = deal(size(V1,1),size(V2,1));
            if sz(1) > sz(2)
                V2 = repelem(V2,sz(1),1);
            end
            for i = 1:sz(1)
                patch(V1{i}(:,1),V1{i}(:,2),V2(i,:),'FaceAlpha',fa,'Linestyle','None');
            end
        end
        %  > 3.2.5. -------------------------------------------------------
        function [] = Var_2D_5(V1,V2,length,c)
            hold on;
            for i = 1:size(V1,1)
                for j = 1:size(V1{i},1)
                    quiver(V1{i}(j,1),V1{i}(j,2),V2{i}(j,1)./length,V2{i}(j,2)./length,'AutoScale','off','Color',c);
                end
            end
        end
        % >> 3.3. ---------------------------------------------------------
        %  > 3.3.1. -------------------------------------------------------
        function [f] = Set_f(msh,ft)
            switch ft
                case "bnd"
                    v = find(~msh.f.logical);
                case "blk"
                    v = find( msh.f.logical);
                otherwise
                    return;
            end
            f = v(randperm(numel(v),1));
        end
        %  > 3.3.2. -------------------------------------------------------
        function [] = SubPlot_pdf(dir,fid)
            exportgraphics(gcf,fid,'ContentType','vector');
            movefile(fid,dir);
        end
    end
end