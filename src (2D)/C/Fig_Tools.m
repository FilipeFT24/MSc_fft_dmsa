classdef Fig_Tools
    methods (Static)         
        %% > Tools.
        % >> 1. -----------------------------------------------------------
        function [] = ChangeLook_1(inp,len,Flag)
            % >> Local variables.
            Xv_i  = inp.msh.lim.Xv_i;
            Xv_f  = inp.msh.lim.Xv_f;
            Yv_i  = inp.msh.lim.Yv_i;
            Yv_f  = inp.msh.lim.Yv_f;
            x_Min = Xv_i-len;
            x_Max = Xv_f+len;
            y_Min = Yv_i-len;
            y_Max = Yv_f+len; 
            
            %  > Other parameters.
            box on; axis equal;
            set(gcf,'color','w');
            set(gca,'Clipping','on');
            set(gca,'Layer','bottom');
            set(gca,'XAxisLocation','bottom');
            set(gca,'YAxisLocation','left');
            set(gca,'TickLabelInterpreter','latex'); set(gca,'FontSize',12);
                       
            %  > X-Axis,Y-Axis label.
            xlabel('$x$','FontSize',20,'Interpreter','latex'); 
            ylabel('$y$','FontSize',20,'Interpreter','latex');
            set(gca,'XLim',[x_Min,x_Max],'XTick',x_Min:(Xv_f-Xv_i)./10:x_Max);
            set(gca,'YLim',[y_Min,y_Max],'YTick',y_Min:(Yv_f-Yv_i)./10:y_Max);
            if Flag
                set(gca,'XTick',[],'YTick',[]);
            end
        end
        % >> 2. -----------------------------------------------------------
        function [] = ChangeLook_2(k,H)
            %  > Other parameters.
            box on; grid on; 
            set(gcf,'color','w');
            set(gca,'GridColor','k');
            set(gca,'GridAlpha',0.25);
            set(gca,'Clipping','on');
            set(gca,'Layer','bottom');
            set(gca,'XAxisLocation','bottom');
            set(gca,'YAxisLocation','left');
            set(gca,'Xdir','reverse'); set(gca,'XScale','log'); set(gca,'YScale','log'); 
            set(gca,'TickLabelInterpreter','latex'); 
            set(gca,'FontSize',10);
            
            %  > X.
            Flag = 'F';
            for i = 1:size(H,2)
                for j = 1:length(H{i})
                    H_Mat(i,j) = H{i}(j);
                end
                if any(size(H{i},2) < 2)
                    Flag = 'T';
                end
            end
            xlabel('$\mathrm{H_{ref}}$','FontSize',12,'Interpreter','latex');
            if strcmpi(Flag,'F')
                xlim([min(H_Mat(:,end)),max(H_Mat(:,1))]);
            end
            %  > Y.
            if k == 1
                ylabel('$|\!|\epsilon_{_{1}}|\!|$','FontSize',10,'Interpreter','latex');
            elseif k == 2
                ylabel('$|\!|\epsilon_{_{2}}|\!|$','FontSize',10,'Interpreter','latex');
            elseif k == 3
                ylabel('$|\!|\epsilon_{_{\infty}}|\!|$','FontSize',10,'Interpreter','latex');
            end
        end       
        % >> 3. -----------------------------------------------------------
        function [ToPatch] = ToPatch_Cell_Face(msh,wdt)
            %  > XYv(1).
            for i = 1:msh.NC+1
                XYv_1(1,i) =  msh.Xv(i);
                XYv_1(2,i) =  wdt;
                XYv_1(3,i) = -wdt;
            end
            %  > XYv(2).           
            XYv_2(1,1)        = 2.*msh.Xv(1)       -msh.Xc(1);
            XYv_2(1,msh.NC+2) = 2.*msh.Xv(msh.NC+1)-msh.Xc(msh.NC);
            for i = 2:msh.NC+1
                XYv_2(2,i) =  wdt;
                XYv_2(3,i) = -wdt;
                if i == 1 || i == msh.NC+2
                    continue;
                else
                    XYv_2(1,i) = msh.Xc(i-1);
                end
            end
            %  > Cell patch.
            for i = 1:msh.NC
                ToPatch.Cell{i}(1,1) = XYv_1(1,i+1);
                ToPatch.Cell{i}(2,1) = XYv_1(2,i+1);
                ToPatch.Cell{i}(1,2) = XYv_1(1,i+1);
                ToPatch.Cell{i}(2,2) = XYv_1(3,i+1);
                ToPatch.Cell{i}(1,3) = XYv_1(1,i);
                ToPatch.Cell{i}(2,3) = XYv_1(3,i);
                ToPatch.Cell{i}(1,4) = XYv_1(1,i);
                ToPatch.Cell{i}(2,4) = XYv_1(2,i);
            end
            %  > Face patch.
            for i = 1:msh.NC+1
                ToPatch.Face{i}(1,1) = XYv_2(1,i+1);
                ToPatch.Face{i}(2,1) = XYv_2(2,i+1);
                ToPatch.Face{i}(1,2) = XYv_2(1,i+1);
                ToPatch.Face{i}(2,2) = XYv_2(3,i+1);
                ToPatch.Face{i}(1,3) = XYv_2(1,i);
                ToPatch.Face{i}(2,3) = XYv_2(3,i);
                ToPatch.Face{i}(1,4) = XYv_2(1,i);
                ToPatch.Face{i}(2,4) = XYv_2(2,i);
            end
        end
        % >> 4. -----------------------------------------------------------
        function [c] = Colormap_cmocean(str)           
            c                      = colorbar();
            c.Location             = 'Eastoutside';
            c.Label.Interpreter    = 'latex';
            c.TickLabelInterpreter = 'latex';
            c.FontSize             =  10;
            c.AxisLocation         = 'out';
            AdvancedColormap(str);
            %  > Colorbar format.
            %  set(c,'xticklabel',cellfun(@(x)sprintf('%.3f',x),num2cell(get(c,'xtick')),'Un',0))
        end
        % >> 5. -----------------------------------------------------------
        function [Marker] = Set_Markers()
            Marker{1} = '-s';
            Marker{2} = '-o';
            Marker{3} = '-^';
            Marker{4} = '-v';
            Marker{5} = '->';
        end
        % >> 7. -----------------------------------------------------------
        function [X_o,Y_o] = Order_Clockwise(Flag,Cell)
            %  > Cell to matrix.
            if strcmpi(Flag,'F')
                Mat = Cell;
            elseif strcmpi(Flag,'T')
                Mat = unique(cell2mat(reshape(Cell,[size(Cell,2),1])),'rows','stable'); 
            end
            X_i = Mat(:,1);
            Y_i = Mat(:,2);
            %  > Unweighted mean.
            Cx = mean(X_i);
            Cy = mean(Y_i);
            %  > Angles.
            th = atan2(Y_i-Cy,X_i-Cx);
            %  > Correct sorted order.
            [~,Order] = sort(th);
            %  > Reorder coordinates:
            X_o = X_i(Order);
            Y_o = Y_i(Order);
        end
        %
        function [] = Plot_Limits(inp,msh,iF)
            % >> Local variables.
            Xv_i = inp.msh.lim.Xv_i;
            Xv_f = inp.msh.lim.Xv_f;
            Yv_i = inp.msh.lim.Yv_i;
            Yv_f = inp.msh.lim.Yv_f;
            
            %  > x_min.
            if round(msh.s.par.lx(1,iF),2) ~= Xv_i
                xline(msh.s.par.lx(1,iF),'-.r','Linewidth',0.5);
            end
            %  > x_max.
            if round(msh.s.par.lx(2,iF),2) ~= Xv_f
                xline(msh.s.par.lx(2,iF),'-.r','Linewidth',0.5);
            end
            %  > y_min.
            if round(msh.s.par.ly(1,iF),2) ~= Yv_i
                yline(msh.s.par.ly(1,iF),'-.r','Linewidth',0.5);
            end
            %  > y_max.
            if round(msh.s.par.ly(2,iF),2) ~= Yv_f
                yline(msh.s.par.ly(2,iF),'-.r','Linewidth',0.5);
            end
        end
        %
        function [] = Export_PDF(Filename,Directory)
            set     (gcf,'PaperSize',[29.7,21.0],'PaperPosition',[0,0,29.7,21.0]);
            print   (gcf,strcat(Filename,'.pdf'),'-dpdf','-r375');
            movefile(strcat(Filename,'.pdf'),Directory); 
        end

    end
end     