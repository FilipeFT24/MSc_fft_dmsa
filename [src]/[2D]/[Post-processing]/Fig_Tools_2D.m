classdef Fig_Tools_2D
    methods (Static)
        %% >> 1. ----------------------------------------------------------
        %  >> 1.1. --------------------------------------------------------
        function [fig] = Set_fig(run,exp)
            fig.run = run;                                     %  > Test working directory.
            fig.exp = exp;                                     %  > Export figure(?).
            fig.C   = linspecer(9,'qualitative');              %  > Colors.
            if ~exp
                fig.LW       =  0.10;                          %  > Line.
                fig.MS       =  7.50;                          %  > Marker size.
                fig.FA       =  0.30;                          %  > FaceAlpha.
                fig.FT_1     = 15.00;                          %  > Legend.
                fig.FT_2     = 14.00;                          %  > x/y-axis.
                fig.FT_3     = [25.00,25.00];                  %  > x/y-label.
                fig.FT_4     = 12.00;                          %  > x/y-axis  (zoom).
                fig.Position = [150,100,1250,600];             %  > Position.
            else
                fig.LW       =  4.00;                          %  > Line.
                fig.MS       =  7.50;                          %  > Marker.
                fig.FA       =  0.25;                          %  > FaceAlpha.
                fig.FT_1     = 22.50;                          %  > Legend.
                fig.FT_2     = 30.00;                          %  > x/y-axis.
                fig.FT_3     = [45.00,45.00];                  %  > x/y-label.
                fig.FT_4     = 20.00;                          %  > x/y-axis  (zoom).
                fig.Position = [350,100,850,600];              %  > Position.
            end
            fig.NT           = [10,10];                        %  > Number of ticks (x/y-direction).
            fig.Folder       = "../[Figures]/[2D]";            %  > Destination folder.
            if ~run
                fig.L{1}     = "$x$";                          %  > x-axis label.
                fig.L{2}     = "$y$";                          %  > y-axis label.
            else
                fig.L{1}     = "$\textrm{NNZ}$";               %  > x-axis label.
                fig.L{2}     = "$\textrm{Error magnitude}$";   %  > y-axis label.
            end
            %  > Markers/line styles.
            %  M = ['+','o','*','x','v','d','^','s','>','<'];
            %  L = ['-','--',':','-.'];
        end
        %  >> 1.2. --------------------------------------------------------
        function [] = ChangeLook(fig,XY)
            %  > Other parameters. 
            axis equal;
            set(gcf,'color','w');
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize',fig.FT_2);
            
            %  > Axis.
            %  > X.
            [x(1),x(2)] = MinMaxElem(XY(:,1)); dx = diff(x)./fig.NT(1); xl = x(1)-dx:dx:x(2)+dx;
            set(gca,'XLim',x,'XTick',xl);
            xlabel(fig.L{1},'FontSize',fig.FT_3(1),'Interpreter','latex');
            %  > Y.
            [y(1),y(2)] = MinMaxElem(XY(:,2)); dy = diff(y)./fig.NT(2); yl = y(1)-dy:dy:y(2)+dy;
            set(gca,'YLim',y,'YTick',yl);
            ylabel(fig.L{2},'FontSize',fig.FT_3(2),'Interpreter','latex');
        end
        %  >> 1.3. --------------------------------------------------------
        function [] = Change_Colormap(fig,str)
            %  > Style.
            ax                     = gca;
            ax_p                   = ax.Position;
            of_p                   = 0.04;
            c                      = colorbar(ax,'Location','east');
            c.Label.Interpreter    = 'latex';
            c.TickLabelInterpreter = 'latex';
            c.FontSize             = fig.FT_2./1.15;
            AdvancedColormap(str);
            %  > Position.
            set(c ,'position',[ax_p(1)+ax_p(3)+of_p,0.165,of_p./2.5,0.705]);
            set(ax,'colorscale','log');
        end
               
        %% >> 2. ----------------------------------------------------------
        %  >> 2.1. --------------------------------------------------------
        function [] = Var_1(f,V,C,M,LW)
            if f
                hold on;
                for i = 1:size(V,1)
                    plot(V{i}(:,1),V{i}(:,2),M,'Color',C,'LineWidth',LW);
                end
            end
        end
        %  >> 2.2. --------------------------------------------------------
        function [] = Var_2(f,V,C,M,MS)
            if f
                hold on;
                for i = 1:size(V,1)
                    plot(V(i,1),V(i,2),M,'Color',C,'MarkerFaceColor',C,'MarkerSize',MS);
                end
            end
        end
        %  >> 2.3. --------------------------------------------------------
        function [] = Var_3(f,V1,V2,FA)
            if f
                hold on;
                [sz(1),sz(2)] = deal(size(V1,1),size(V2,1));
                if sz(1) > sz(2)
                    V2 = repelem(V2,sz(1),1);
                end
                for i = 1:sz(1)
                    patch(V1{i}(:,1),V1{i}(:,2),V2(i,:),'FaceAlpha',FA,'Linestyle','None');
                end
            end
        end
        %  >> 2.4. --------------------------------------------------------
        function [] = Var_4(f,V1,V2,x,C)
            if f
                hold on;
                for i = 1:size(V1,1)
                    for j = 1:size(V1{i},1)
                        quiver(V1{i}(j,1),V1{i}(j,2),V2{i}(j,1)./x,V2{i}(j,2)./x,'AutoScale','off','Color',C);
                    end
                end
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
                Fig_Tools_2D.Export_PDF(fid,fig.Folder);
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