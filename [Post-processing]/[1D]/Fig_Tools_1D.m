classdef Fig_Tools_1D
    methods (Static)
        %% >> 1. ----------------------------------------------------------
        %  >> 1.1. --------------------------------------------------------
        function [fig] = Set_fig(run,exp)
            fig.run = run;                                    %  > Test working directory.
            fig.exp = exp;                                    %  > Export figure(?).
            if ~exp
                fig.LW       =  2.00;                         %  > Line.
                fig.MS       =  3.50;                         %  > Marker size.
                fig.FT_1     = 15.00;                         %  > Legend.
                fig.FT_2     = 15.00;                         %  > x/y-axis.
                if ~run
                    fig.FT_3 = [22.50,17.50];                 %  > x/y-label.
                else
                    fig.FT_3 = [15.00,17.50];                 %  > x/y-label.
                end
                fig.FT_4     = 11.00;                         %  > x/y-axis  (zoom).
                fig.Position = [150,100,1250,600];            %  > Position.
            else
                fig.LW       =  4.00;                         %  > Line.
                fig.MS       =  7.50;                         %  > Marker.
                fig.FT_1     = 22.50;                         %  > Legend.
                fig.FT_2     = 30.00;                         %  > x/y-axis.
                if ~run
                    fig.FT_3     = [45.00,35.00];             %  > x/y-label.
                else
                    fig.FT_3     = [35.00,35.00];             %  > x/y-label.
                end
                fig.FT_4     = 20.00;                         %  > x/y-axis  (zoom).
                fig.Position = [350,100,850,600];             %  > Position.
            end
            fig.trsh     = 1.0e-12;                           %  > Do not plot below 'trsh'.
            fig.nsh      = 0;                                 %  > Number of elements below 'trsh'.
            fig.NT       = [10,10];                           %  > Number of ticks (x/y-direction).
            fig.Folder   = "../[Figures]/[1D]";               %  > Destination folder.
            if ~run
                fig.c       = 5.0e-03;                        %  > x-axis width.
                fig.L{1}    = "$x$";                          %  > x-axis label.
                fig.L{2}    = "$\textrm{Error magnitude}$";   %  > y-axis label.
            else
                fig.c       = 0;                              %  > x-axis width.
                fig.L{1}    = "$\textrm{NNZ}$";               %  > x-axis label.
                fig.L{2}    = "$\textrm{Error magnitude}$";   %  > y-axis label.
            end
            %  > Markers/line styles.
            %  M = ['+','o','*','x','v','d','^','s','>','<'];
            %  L = ['-','--',':','-.'];
        end
        %  >> 1.2. --------------------------------------------------------
        function [] = ChangeLook_1D(fig,XM,L_XY)
            %  > Other parameters.
            [xl(1),xl(2)] = MinMaxElem(XM);
            box on;
            set(gcf,'color','w');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fig.FT_2);
            
            % >> Axis.
            if ~fig.run
                %  > X-Axis.
                for i = 1:fig.NT(1)+1
                    xt(i) = xl(1)-(xl(1)-xl(2))./fig.NT(1).*(i-1);
                end
                xticks(xt); 
            end
            xlim([xl(1)-fig.c,xl(2)+fig.c]);
            xlabel(L_XY(1),'FontSize',fig.FT_3(1),'Interpreter','latex');
            %  > Y-Axis.
            ylabel(L_XY(2),'FontSize',fig.FT_3(2),'Interpreter','latex');
        end
        %% >> 2. ----------------------------------------------------------
        %  >> 2.1. --------------------------------------------------------
        function [LY,P,YV] = Var_1(fig,M,LX,X,Y)
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
                    P {k}             = plot(X(:,i),Y(:,i),M(i),'Color',fig.C(i,:),'LineWidth',fig.LW,'MarkerFaceColor',fig.C(i,:),'MarkerSize',fig.MS);
                    [YV(k,1),YV(k,2)] = MinMaxElem(Y(j{i},i));
                end
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [LY,P,YV] = Var_2(fig,M,LX,X,Y)
            hold on;
            k = 0;
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
        % >> 2.3. ---------------------------------------------------------
        %  > 2.3.1. -------------------------------------------------------
        function [] = Set_Plot_1(fig,y,XM)
            Fig_Tools_1D.ChangeLook_1D(fig,XM,fig.L);
            ax = gca; box off;
            dX = [XM(1),XM(end)];
            dY = [-y,y];
            xlim(dX); xticks(dX); set(gca,'XTickLabel',[]);
            ylim(dY); yticks(dY); set(gca,'YTickLabel',[]); ax.YAxis.Visible = 'off';
        end
        %  > 2.3.2. -------------------------------------------------------
        function [] = Set_Plot_2(fig,L,P,XM,YM_1,YM_2,NC)
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',fig.FT_1,'NumColumns',NC);
            if fig.run
                set(gca,'XScale','log');
            end
            set(gca,'YScale','log');
            ylim([10.^(ceil(log10(min(YM_1(:,1))))+YM_2(1)),...
                  10.^(ceil(log10(max(YM_1(:,2))))+YM_2(2))]);
            Fig_Tools_1D.ChangeLook_1D(fig,XM,fig.L);
        end
        %  > 2.3.3. -------------------------------------------------------
        function [str] = Set_str_1(i)
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
        % >> 3. -----------------------------------------------------------
        function [] = Exp_Or_Zoom(fig,zoom,fid)
            %  > Zoom(?).
            if zoom
                zp = BaseZoom;
                zp.plot(fig);
            end
            %  > Export(?).
            if fig.exp
                Fig_Tools_1D.Export_PDF(fid,fig.Folder);
            end
        end
        % >> 4. -----------------------------------------------------------
        function [] = Export_PDF(Filename,Directory)
            set     (gcf,'PaperSize',[29.7,21.5],'PaperPosition',[0,0,29.7,21.5]);
            print   (gcf,strcat(Filename,'.pdf'),'-dpdf','-r500');
            movefile(    strcat(Filename,'.pdf'),Directory);
        end
    end
end