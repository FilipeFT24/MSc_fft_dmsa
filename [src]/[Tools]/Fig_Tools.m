classdef Fig_Tools
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [fig] = Set_fig(exp,run,zoom)
            fig.exp  = exp;                                  %  > Export figure(?).
            fig.run  = run;                                  %  > Test working directory.
            fig.zoom = zoom;                                 %  > Zoom(?).
            fig.C    = linspecer(9,'qualitative');           %  > Colors.
            if ~exp
                fig.LW       =  1.50;                        %  > Line.
                fig.MS       =  5.00;                        %  > Marker size.
                fig.FA       =  0.25;                        %  > FaceAlpha.
                fig.FT_1     = 15.00;                        %  > Legend.
                fig.FT_2     = 14.00;                        %  > x/y-axis.
                fig.FT_3     = [25.00,25.00];                %  > x/y-label.
                fig.FT_4     = 12.50;                        %  > x/y-axis (zoom).
                fig.Position = [150,100,1250,600];           %  > Position.
            else
                fig.LW       =  0.10;                        %  > Line.
                fig.MS       =  7.50;                        %  > Marker.
                fig.FA       =  0.25;                        %  > FaceAlpha.
                fig.FT_1     = 25.00;                        %  > Legend.
                fig.FT_2     = 30.00;                        %  > x/y-axis.
                fig.FT_3     = [42.50,42.50];                %  > x/y-label.
                fig.FT_4     = 25.00;                        %  > x/y-axis (zoom).
                fig.Position = [350,100,850,600];            %  > Position.
            end
            fig.nsh          = 0;                            %  > Number of elements below "trsh".
            fig.trsh         = 1.0e-16;                      %  > Error treshold.
            fig.NT           = [10,10];                      %  > Number of ticks (x/y-direction).
            fig.Folder       = "../[Figures]";               %  > Destination folder.
            if ~run
                fig.L{1}     = "$x$";                        %  > x-axis label.
                fig.L{2}     = "$y$";                        %  > y-axis label.
            else
                fig.L{1}     = "$\textrm{NNZ}$";             %  > x-axis label.
                fig.L{2}     = "$\textrm{Error magnitude}$"; %  > y-axis label.
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
                    str = "\phantom{\nabla}\phi";
                case 2
                    str = "\nabla\phi";
                otherwise
                    return;
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        function [str] = Set_str_2(i)
            switch i
                case 1
                    str = "\phi";
                case 2
                    str = "\nabla\phi";
                otherwise
                    return;
            end
        end
        %  > 1.2.3. -------------------------------------------------------
        function [str] = Set_str_3(i)
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
        function [] = Map_1D_1(fig,fx,XM,L_XY)
            %  > Other parameters.
            [xl(1),xl(2)] = MinMaxElem(XM);
            box on;
            set(gcf,'color','w');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fig.FT_2);
            
            %  > Axis.
            if fx
                i     = 1:fig.NT(1)+1;
                xt(i) = xl(1)-(xl(1)-xl(2))./fig.NT(1).*(i-1);
                xticks(xt);
            end
            xlim  ([xl(1),xl(2)]);
            xlabel(L_XY(1),'FontSize',fig.FT_3(1),'Interpreter','latex');
            ylabel(L_XY(2),'FontSize',fig.FT_3(2),'Interpreter','latex');
        end
        %  > 2.1.2. -------------------------------------------------------
        function [] = Map_1D_2(fig,L,P,fx,XM,YM_1,YM_2,NC)
            ax = gca;
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',fig.FT_1,'NumColumns',NC);
            if fig.run
                set(ax,'XScale','log'); ax.XAxis.Exponent = 3;
            end
            set(ax,'YScale','log');
            if ~isempty(YM_1) && ~isempty(YM_2)
                a = 10.^(ceil(log10(min(YM_1(:,1))))+YM_2(1));
                b = 10.^(ceil(log10(max(YM_1(:,2))))+YM_2(2));
                if a < fig.trsh
                    a = fig.trsh;
                end
                ylim([a,b]);
            end
            Fig_Tools.Map_1D_1(fig,fx,XM,fig.L);
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        function [LY,P,YV] = Var_1D_1(fig,M,LX,X,Y)
            %  > Auxiliary variables.
            a = size(X,2);
            b = size(Y,2);
            if a < b
                X = repelem(X,1,b);
            end
            
            hold on;
            k = 0;
            for i = 1:b
                j{i} = Y(:,i) > fig.trsh; Y(~j{i},i) = 0;
                if nnz(j{i})  > fig.nsh
                    k                 = k+1;
                    LY{k}             = LX{i};
                    P {k}             = plot(X(:,i),Y(:,i),M(i),'Color',fig.C(k,:),'LineWidth',fig.LW,'MarkerFaceColor',fig.C(k,:),'MarkerSize',fig.MS);
                    [YV(k,1),YV(k,2)] = MinMaxElem(Y(j{i},i));
                end
            end
        end
        
        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > 3.1.1. -------------------------------------------------------
        function [] = Map_2D_1(msh,fig)
            %  > Other parameters.
            axis equal;
            set(gcf,'color','w');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fig.FT_2);
            
            %  > Axis.
            XY = cat(1,msh.f.xy.v{:});
            %  > X.
            [x(1),x(2)] = MinMaxElem(XY(:,1)); dx = diff(x)./fig.NT(1); xl = x(1)-dx:dx:x(2)+dx;
            set(gca,'XLim',x,'XTick',xl);
            xlabel(fig.L{1},'FontSize',fig.FT_3(1),'Interpreter','latex');
            %  > Y.
            [y(1),y(2)] = MinMaxElem(XY(:,2)); dy = diff(y)./fig.NT(2); yl = y(1)-dy:dy:y(2)+dy;
            set(gca,'YLim',y,'YTick',yl);
            ylabel(fig.L{2},'FontSize',fig.FT_3(2),'Interpreter','latex');
            %  > Grid.
            LW = 0.1;
            Fig_Tools.Var_2D_1(msh.f.xy.v,"k","-",LW);
        end
        %  > 3.1.2. -------------------------------------------------------
        function [] = Map_2D_2(fig,y)
            %  > Style.
            a                      = gca;
            p                      = a.Position;
            of_p                   = 0.035;
            c                      = colorbar(a,'Location','east');
            c.Label.Interpreter    = 'latex';
            c.TickLabelInterpreter = 'latex';
            c.FontSize             = fig.FT_1;
            AdvancedColormap(y.str);
            %  > Limits.
            if y.c(1) >= fig.trsh
                y_lim_1 = 10.^(ceil(log10(y.c(1)))-1);
            else
                y_lim_1 = fig.trsh;
            end
            caxis([y_lim_1,10.^(ceil(log10(y.c(2)))+0)]);
            %  > Position.
            set(a,'colorscale','log');
            set(c,'position',[p(1)+p(3)+of_p,0.165,of_p./2.5,0.705]);
            %  > Title.
            title(y.title,'interpreter','latex');
        end
        %  > 3.1.3. -------------------------------------------------------
        function [] = Map_2D_3(fig)
            if fig.zoom
                zp = BaseZoom;
                zp.plot(fig);
            end
            if fig.exp
                set     (gcf,'PaperSize',[29.7,21.5],'PaperPosition',[0,0,29.7,21.5]);
                print   (gcf,strcat(Filename,'.pdf'),'-dpdf','-r500');
                movefile(    strcat(Filename,'.pdf'),Directory);
            end
        end
        % >> 3.2. ---------------------------------------------------------
        %  > 3.2.1. -------------------------------------------------------
        function [] = Var_2D_1(V,C,M,LW)
            hold on;
            for i = 1:size(V,1)
                plot(V{i}(:,1),V{i}(:,2),M,'Color',C,'LineWidth',LW);
            end
        end
        %  > 3.2.2. -------------------------------------------------------
        function [P] = Var_2D_2(V,C,M,MS)
            hold on;
            P = plot(V(:,1),V(:,2),M,'Color',C,'MarkerFaceColor',C,'MarkerSize',MS);
        end
        %  > 3.2.3. -------------------------------------------------------
        function [] = Var_2D_3(V1,V2,FA)
            hold on;
            [sz(1),sz(2)] = deal(size(V1,1),size(V2,1));
            if sz(1) > sz(2)
                V2 = repelem(V2,sz(1),1);
            end
            for i = 1:sz(1)
                patch(V1{i}(:,1),V1{i}(:,2),V2(i,:),'FaceAlpha',FA,'Linestyle','None');
            end
        end
        %  > 3.2.4. -------------------------------------------------------
        function [] = Var_2D_4(V1,V2,x,C)
            hold on;
            for i = 1:size(V1,1)
                for j = 1:size(V1{i},1)
                    quiver(V1{i}(j,1),V1{i}(j,2),V2{i}(j,1)./x,V2{i}(j,2)./x,'AutoScale','off','Color',C);
                end
            end
        end
    end
end