classdef Fig_Tools_1D
    methods (Static)         
        %% > Tools.
        % >> 1. -----------------------------------------------------------
        function [] = ChangeLook_1D(box_w,xc,NX,sz)
            % >> Other parameters.
            if box_w
                box on;
            end
            set(gcf,'color','w');
            set(gca,'Clipping','on');
            set(gca,'Layer','bottom');
            set(gca,'XAxisLocation','bottom');
            set(gca,'TickLength',[0,0]);
            set(gca,'XTick',[]);
            set(gca,'TickLabelInterpreter','latex'); 
            set(gca,'FontSize',sz);
            % >> Labels.
            %  > X-Axis label.
            xlim_label = [min(xc),max(xc)];
            for ix_ticks = 1:NX+1
                xticks_label(ix_ticks) = xlim_label(1)-(xlim_label(1)-xlim_label(2))./NX.*(ix_ticks-1);
            end
            xlabel('$x$','FontSize',20,'Interpreter','latex'); xlim([xlim_label(1),xlim_label(2)]); xticks(xticks_label);
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
        function [c] = Colormap_style(str_1,str_2,yv,NY,sz) %#ok<INUSL>
            c                      = colorbar();
            c.Location             = 'Eastoutside';
            c.Label.Interpreter    = 'latex';
            c.TickLabelInterpreter = 'latex';
            c.FontSize             =  sz;
            c.AxisLocation         = 'out';
            c.Label.String         = str_1;
            AdvancedColormap(str_2);
            %  > Colorbar format.
            %  c.Ticks = linspace(yv(1),yv(2),NY);
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