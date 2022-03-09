classdef Fig_Tools_1D
    methods (Static)         
        %% > Tools.
        % >> 1. -----------------------------------------------------------
        function [] = ChangeLook_1D(box_w,x_w,y_w,xc,NX,L_X,L_Y,SZ_X,SZ_Y)
            %  > Other parameters.
            if box_w
                box on;
            end
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
                c = 0.0025;
                xlabel(L_X,'FontSize',SZ_X.*1.50,'Interpreter','latex'); xlim([xlim_label(1)-c,xlim_label(2)+c]); xticks(xticks_label);
            end
            %  > Y-Axis.
            if y_w
                ylabel(L_Y,'FontSize',SZ_Y.*1.25,'Interpreter','latex');
            end
        end
        % >> 2. -----------------------------------------------------------
        function [fig] = Set_fig(Exp)
            if ~Exp
                fig.LW_1 =  1.5;
                fig.LW_2 =  1.5;
                fig.MS_1 =  3.0;
                fig.FT_1 = 14.0;
                fig.FT_2 = 12.5;
                fig.FT_3 = 12.5;
            else
                fig.LW_1 =  3.0;
                fig.LW_2 =  2.5;
                fig.MS_1 =  5.0;
                fig.FT_1 = 25.0;
                fig.FT_2 = 30.0;
                fig.FT_3 = 30.0;
            end
        end
        % >> 3. -----------------------------------------------------------
        function [] = Export_PDF(Filename,Directory)
            set     (gcf,'PaperSize',[28.50,20.85],'PaperPosition',[0,0,29.7,21.5]);
            print   (gcf,strcat(Filename,'.pdf'),'-dpdf','-r500');
            movefile(strcat(Filename,'.pdf'),Directory); 
        end
    end
end   