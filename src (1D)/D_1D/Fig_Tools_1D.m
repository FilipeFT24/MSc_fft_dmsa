classdef Fig_Tools_1D
    methods (Static)
        %% >> 1. ----------------------------------------------------------
        %  >> 1.1. --------------------------------------------------------
        function [fig] = Set_fig(Exp)
            if ~Exp
                fig.LW       =  1.25;               %  > Line.
                fig.MS       =  4.00;               %  > Marker size.
                fig.FT_1     = 14.00;               %  > Legend.
                fig.FT_2     = 15.00;               %  > x/y-axis.
                fig.FT_3     = [22.50,17.50];       %  > x/y-label.
                fig.FT_4     = 11.00;               %  > x/y-axis  (zoom).
                fig.Position = [150,100,1250,600];  %  > Position.
            else
                fig.LW       =  3.50;               %  > Line.
                fig.MS       =  5.00;               %  > Marker.
                fig.FT_1     = 25.00;               %  > Legend.
                fig.FT_2     = 30.00;               %  > x/y-axis.
                fig.FT_3     = [42.50,35.00];       %  > x/y-label.
                fig.FT_4     = 20.00;               %  > x/y-axis  (zoom).
                fig.Position = [350,100,850,600];   %  > Position.  
            end
            fig.c        = 2.5e-03;                 %  > x-axis width.
            fig.trsh     = 1.0e-10;                 %  > Do not plot below 'trsh'.
            fig.nsh      = 5;                       %  > Number of elements below 'trsh'.
            fig.NT       = [10,10];                 %  > Number of ticks (x/y-direction).
            fig.Folder   = "../[Figures]/[1D]";     %  > Destination folder.
            fig.Label{1} = "$x$";
            fig.Label{2} = "$\textrm{Error magnitude}$";
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
        function [LY,P,YV] = Var_1(fig,LX,X,Y)
            hold on;
            k = 0;
            for i = 1:size(Y,2)
                j{i} = Y(:,i) > fig.trsh; Y(~j{i},i) = 0;
                if nnz(j{i})  > fig.nsh
                    k                 = k+1;
                    LY{k}             = LX{i};
                    P {k}             = plot(X,Y(:,i),fig.M(i),'Color',fig.C(i,:),'LineWidth',fig.LW,'MarkerFaceColor',fig.C(i,:),'MarkerSize',fig.MS);
                    [YV(k,1),YV(k,2)] = MinMaxElem(Y(j{i},i));
                end
            end
        end
        % >> 2.2. ---------------------------------------------------------
        function [LY,P,YV] = Var_2(fig,LX,X,Y)
            hold on;
            k = 0;
            m = 1;
            for i = 1:size(Y,2)
                j = Y(m,i) > fig.trsh; Y(~j,i) = 0;
                if j
                    k                 = k+1;
                    P {k}             = line([X(1),X(end)],[Y(m,i),Y(m,i)],'Color',fig.C(i,:),'Linewidth',fig.LW,'Linestyle','-');
                    LY{k}             = LX{i};
                    [YV(k,1),YV(k,2)] = MinMaxElem(Y(m,i));
                end
            end
        end
        % >> 2.3. ---------------------------------------------------------
        %  > 2.3.1. -------------------------------------------------------
        function [] = Set_Plot(fig,msh,L,P,YM,NC)
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',fig.FT_1,'NumColumns',NC);
            set(gca,'YScale','log');
            ylim([10.^(ceil(log10(min(YM(:,1))))-1),...
                  10.^(ceil(log10(max(YM(:,2))))+1)]);
            Fig_Tools_1D.ChangeLook_1D(fig,msh,fig.Label);
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
        function [] = Exp_Or_Zoom(edt,fid,folder)
            %  > Export(?).
            if edt(1)
                Fig_Tools_1D.Export_PDF(fid,folder)
            end
            %  > Zoom(?).
            if edt(2)
                zp = BaseZoom;
                zp.plot(edt(1));
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