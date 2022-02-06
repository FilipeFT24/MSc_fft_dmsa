classdef Fig_Tools_1D
    methods (Static)         
        %% > Tools.
        % >> 1. -----------------------------------------------------------
        function [] = ChangeLook_1D(box_w,y_w,xc,NX,L_X,L_Y,SZ_X,SZ_Y)
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
            xlim_label = [min(xc),max(xc)];
            for ix_ticks = 1:NX+1
                xticks_label(ix_ticks) = xlim_label(1)-(xlim_label(1)-xlim_label(2))./NX.*(ix_ticks-1);
            end
            set   (gca,'XAxisLocation','bottom');
            set   (gca,'XTick',[]);
            xlabel(L_X,'FontSize',SZ_X,'Interpreter','latex'); xlim([xlim_label(1),xlim_label(2)]); xticks(xticks_label);
            %  > Y-Axis.
            if y_w
                ylabel(L_Y,'FontSize',SZ_Y,'Interpreter','latex');
            end
        end
        % >> 2. -----------------------------------------------------------
        % >> 3. -----------------------------------------------------------
        function [c_xy] = ToPatch(msh,wdt)
            %  > XY_v.
            i       =  1:msh.f.NF;
            v_xy(i) =  msh.f.Xv(i);
            %  > c_xy.
            for i = 1:msh.c.NC
                c_xy{i}(1,[1,2]) =  v_xy(i+1);
                c_xy{i}(1,[3,4]) =  v_xy(i);
                c_xy{i}(2,[1,4]) =  wdt;
                c_xy{i}(2,[2,3]) = -wdt;
            end
        end
        % >> 4. -----------------------------------------------------------
        function [c] = Colormap_style(L1,L2,SZ_Y)
            c                      =  colorbar();
            c.Location             = 'Eastoutside';
            c.Label.Interpreter    = 'latex';
            c.TickLabelInterpreter = 'latex';
            c.FontSize             =  SZ_Y;
            c.AxisLocation         = 'out';
            c.Label.String         =  L1;
            AdvancedColormap(L2);
            %  > Colorbar format.
            %  set(c,'xticklabel',cellfun(@(x)sprintf('%.3f',x),num2cell(get(c,'xtick')),'Un',0))
        end
        % >> 5. -----------------------------------------------------------
        % >> 6. -----------------------------------------------------------
        function [] = Export_PDF(Filename,Directory)
            set     (gcf,'PaperSize',[29.7,21.0],'PaperPosition',[0,0,29.7,21.0]);
            print   (gcf,strcat(Filename,'.pdf'),'-dpdf','-r500');
            movefile(strcat(Filename,'.pdf'),Directory); 
        end
    end
end   