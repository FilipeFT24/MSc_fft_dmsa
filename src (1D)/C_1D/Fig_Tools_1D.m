classdef Fig_Tools_1D
    methods (Static)
        %% >> 1. ----------------------------------------------------------
        %  >> 1.1. --------------------------------------------------------
        function [fig] = Set_fig(Exp)
            if ~Exp
                fig.LW       =  2.0;                %  > Line.
                fig.MS       =  3.0;                %  > Marker size.
                fig.FT_1     = 12.5;                %  > Legend.
                fig.FT_2     = 15.0;                %  > x/y-axis.
                fig.FT_3     = [22.5,17.5];         %  > x/y-label.
                fig.Position = [150,100,1250,600];  %  > Position.
            else
                fig.LW       =  3.5;                %  > Line.
                fig.MS       =  5.0;                %  > Marker.
                fig.FT_1     = 25.0;                %  > Legend.
                fig.FT_2     = 30.0;                %  > x/y-axis.
                fig.FT_3     = [42.5,35.0];         %  > x/y-label.
                fig.Position = [350,100,850,600];   %  > Position.
                fig.Folder   = "../[Figures]/[1D]"; %  > Destination folder.
            end
            fig.c    = 2.5e-03;                     %  > x-axis width.
            fig.trsh = 1.0e-12;                     %  > Do not plot below 'trsh'.
            fig.nsh  = 10;                          %  > Number of elements below 'trsh'.
            fig.NT   = [10,10];                     %  > Number of ticks (x/y-direction).
            %  > Markers/line styles.
            %  M = ['+','o','*','x','v','d','^','s','>','<'];
            %  L = ['-','--',':','-.'];
        end
        %  >> 1.2. --------------------------------------------------------
        function [] = ChangeLook_1D(fig,msh,L_XY)
            %  > Other parameters.
            [xl(1),xl(2)] = MinMaxElem(msh.f.Xv);
            box on;
            set(gcf,'color','w');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fig.FT_2);
            
            % >> Axis.
            %  > X-Axis.
            for i = 1:fig.NT(1)+1
                xt(i) = xl(1)-(xl(1)-xl(2))./fig.NT(1).*(i-1);
            end
            xticks(xt);
            xlim([xl(1)-fig.c,xl(2)+fig.c]);
            xlabel(L_XY(1),'FontSize',fig.FT_3(1),'Interpreter','latex');
            %  > Y-Axis.
            ylabel(L_XY(2),'FontSize',fig.FT_3(2),'Interpreter','latex');
        end
        %% >> 2. ----------------------------------------------------------
        %  >> 2.1. --------------------------------------------------------
        function [fig,P,YV] = Var_1(fig,X,Y)
            hold on;
            k = 0;
            for i = 1:size(Y,2)
                j{i} = Y(:,i) > fig.trsh; Y(~j{i},i) = 0;
                if nnz(j{i})  > fig.nsh
                    k                 = k+1;
                    P{k}              = plot(X,Y(:,i),fig.M(i),'Color',fig.C(i,:),'LineWidth',fig.LW,'MarkerFaceColor',fig.C(i,:),'MarkerSize',fig.MS);
                    L{k}              = fig.L1{i};
                    [YV(k,1),YV(k,2)] = MinMaxElem(Y(j{i},i));
                end
            end
            fig.L1 = L;
        end
        % >> 2.2. ---------------------------------------------------------
        function [fig,P,YV] = Var_2(fig,X,Y)
            hold on;
            k = 0;
            m = 1;
            for i = 1:size(Y,2)
                j = Y(m,i) > fig.trsh; Y(~j,i) = 0;
                if j
                    k                 = k+1;
                    P{k}              = line([X(1),X(end)],[Y(m,i),Y(m,i)],'Color',fig.C(i,:),'Linewidth',fig.LW,'Linestyle','-');
                    L{k}              = fig.L1{i};
                    [YV(k,1),YV(k,2)] = MinMaxElem(Y(m,i));
                end
            end
            fig.L1 = L;
        end
        % >> 2.3. ---------------------------------------------------------
        %  > 2.3.1. -------------------------------------------------------
        function [] = Set_Plot(fig,msh,P,YM)
            legend([P{:}],[fig.L1{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',fig.FT_1,'NumColumns',2);
            set(gca,'YScale','log');
            ylim([10.^(ceil(log10(min(YM(:,1))))-1),...
                  10.^(ceil(log10(max(YM(:,2))))+1)]);
            Fig_Tools_1D.ChangeLook_1D(fig,msh,fig.L2);
        end
        %  > 2.3.2. -------------------------------------------------------
        function [str] = Switch_Legend(i)
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
        function [] = Export_PDF(Filename,Directory)
            set     (gcf,'PaperSize',[29.7,21.5],'PaperPosition',[0,0,29.7,21.5]);
            print   (gcf,strcat(Filename,'.pdf'),'-dpdf','-r500');
            movefile(strcat(Filename,'.pdf'),Directory);
        end
    end
end