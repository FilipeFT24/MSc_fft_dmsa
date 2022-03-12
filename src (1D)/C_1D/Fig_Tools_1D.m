classdef Fig_Tools_1D
    methods (Static)         
        %% > Tools.
        % >> 1. -----------------------------------------------------------
        function [] = ChangeLook_1D(x_w,y_w,xc,NX,L_X,L_Y,SZ_X,SZ_Y)
            %  > Other parameters.
            box on;
            set(gcf,'color','w');
            set(gca,'Clipping','on');
            set(gca,'Layer','bottom');
            set(gca,'TickLabelInterpreter','latex'); 
            set(gca,'FontSize',min(SZ_X,SZ_Y));
            %  > X-Axis.
            if x_w
                set(gca,'xtick',[])
                xlim_label = [min(xc),max(xc)];
                for ix_ticks = 1:NX+1
                    xticks_label(ix_ticks) = xlim_label(1)-(xlim_label(1)-xlim_label(2))./NX.*(ix_ticks-1);
                end
                c = 2.5e-3;
                xlabel(L_X,'FontSize',SZ_X,'Interpreter','latex'); xlim([xlim_label(1)-c,xlim_label(2)+c]); xticks(xticks_label);
            end
            %  > Y-Axis.
            if y_w
                ylabel(L_Y,'FontSize',SZ_Y,'Interpreter','latex');
            end
        end
        % >> 2. -----------------------------------------------------------
        function [fig] = Set_fig(Exp)
            if ~Exp
                fig.LW   =  2.0; %  > Line.
                fig.MS   =  3.0; %  > Marker size.
                fig.FT_1 = 13.0; %  > Legend.
                fig.FT_2 = 15.0; %  > x-label/y-label.
            else
                fig.LW   =  3.0; %  > Line.
                fig.MS   =  5.0; %  > Marker.
                fig.FT_1 = 25.0; %  > Legend.
                fig.FT_2 = 30.0; %  > x-label/y-label.
            end
            fig.NT   = 10;       %  > Number of ticks (x-direction).
            fig.MT   = [":o","-.v","-s","--s"];
            fig.trsh = 10e-12;
        end
        % >> 3. -----------------------------------------------------------
        %  > 3.1. ---------------------------------------------------------
        %  > Plot variable #1.
        function [P,YV] = Var_1(fig,C,X,Y)
            hold on;
            for i = 1:size(Y,2)
                j                 = Y(:,i) > fig.trsh;
                P{i}              = plot(X,Y(:,i),':o','Color',C(i,:),'LineWidth',fig.LW,'MarkerFaceColor',C(i,:),'MarkerSize',fig.MS);
                [YV(i,1),YV(i,2)] = MinMaxElem(Y(j,i)); 
            end
        end
        %  > 3.2. ---------------------------------------------------------
        %  > Plot variable #2.
        function [P,YV] = Var_2(fig,C,X,Y)
            hold on;
            for i = 1:size(Y,2)
                j                 = Y(i) > fig.trsh;
                P{i}              = line([X(1),X(end)],[Y(i),Y(i)],'Color',C(i,:),'Linewidth',fig.LW,'Linestyle','-');
                [YV(i,1),YV(i,2)] = MinMaxElem(Y(j,i)); 
            end
        end
        %  > 3.3. ---------------------------------------------------------
        function [] = Set_Plot(fig,msh,L,P,YM)
            legend([P{:}],[L{:}],...
                'Interpreter','latex','Location','Northeast','FontSize',fig.FT_1,'NumColumns',2);
            set(gca,'YScale','log');
            ylim([10.^(ceil(log10(min(YM(:,1))))-1),...
                  10.^(ceil(log10(max(YM(:,2))))+1)]);
            Fig_Tools_1D.ChangeLook_1D(1,1,msh.f.Xv,fig.NT,"$x$",L{end},fig.FT_2,fig.FT_2);
        end
        
        % >> 4. -----------------------------------------------------------
        function [] = Export_PDF(Filename,Directory)
            set     (gcf,'PaperSize',[28.50,20.85],'PaperPosition',[0,0,29.7,21.5]);
            print   (gcf,strcat(Filename,'.pdf'),'-dpdf','-r500');
            movefile(strcat(Filename,'.pdf'),Directory);
        end
    end
end   