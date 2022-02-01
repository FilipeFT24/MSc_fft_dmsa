classdef Fig_Tools_1D
    methods (Static)         
        %% > Tools.
        %  > #1.
        function [] = ChangeLook_1(x_Min,x_Max,NX)
            %  > Other parameters.
            box on; 
            set(gcf,'color','w');
            set(gca,'Clipping','on');
            set(gca,'Layer','bottom');
            set(gca,'XAxisLocation','bottom');
            set(gca,'YAxisLocation','left');
            
            %  > X-Axis,Y-Axis label.
            xlim_label = [x_Min,x_Max];
            for i_ticks = 1:NX+1
                xticks_label(i_ticks) = xlim_label(1)-(xlim_label(1)-xlim_label(2))./NX.*(i_ticks-1);
            end
            set(gca,'TickLabelInterpreter','latex'); set(gca,'FontSize',10);
            xlabel('x','FontSize',10,'Interpreter','latex'); xlim([xlim_label(1),xlim_label(2)]); xticks(xticks_label);
        end
        %  > #2.
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
            
            %  > x.
            Flag = 'F';
            for i = 1:size(H,2)
                for j = 1:length(H{i})
                    H_Mat(i,j) = H{i}(j);
                end
                if any(size(H{i},2) < 2)
                    Flag = 'T';
                end
            end
            xlabel('$\Delta\mathrm{x}$','FontSize',10,'Interpreter','latex');
            if strcmpi(Flag,'F')
                xlim([min(H_Mat(:,end)),max(H_Mat(:,1))]);
            end
            %  > y.
            if k == 1
                ylabel('$|\!|\epsilon_{_{1}}|\!|$','FontSize',10,'Interpreter','latex');
            elseif k == 2
                ylabel('$|\!|\epsilon_{_{2}}|\!|$','FontSize',10,'Interpreter','latex');
            elseif k == 3
                ylabel('$|\!|\epsilon_{_{\infty}}|\!|$','FontSize',10,'Interpreter','latex');
            end
        end       
        % >> #3.
        function [ToPatch] = ToPatch_Cell_Face(msh,wdt)
            % >> XYv(1).
            for i = 1:msh.c.NC+1
                XYv_1(1,i) =  msh.f.Xv(i);
                XYv_1(2,i) =  wdt;
                XYv_1(3,i) = -wdt;
            end
            % >> XYv(2).           
            XYv_2(1,1)        = 2.*msh.f.Xv(1)       -msh.c.Xc(1);
            XYv_2(1,msh.c.NC+2) = 2.*msh.f.Xv(msh.c.NC+1)-msh.c.Xc(msh.c.NC);
            for i = 2:msh.c.NC+1
                XYv_2(2,i) =  wdt;
                XYv_2(3,i) = -wdt;
                if i == 1 || i == msh.c.NC+2
                    continue;
                else
                    XYv_2(1,i) = msh.c.Xc(i-1);
                end
            end
            %  > Cell patch.
            for i = 1:msh.c.NC
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
            for i = 1:msh.c.NC+1
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
        % >> #4.
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
        % >> #5.
        function [Marker] = Set_Markers()
            Marker{1} = '-s';
            Marker{2} = '-o';
            Marker{3} = '-^';
            Marker{4} = '-v';
            Marker{5} = '->';
        end
        % >> #6.
        function [H_av,CR] = Compute_ConvergenceRate(H,X)
            for i = 1:size(X,2)-1
                for j = 1:3
                    CR  (i,j) = log(cell2mat(X{i+1}(j))./cell2mat(X{i}(j)))./log(H(i+1)./H(i));
                    H_av(i)   = 1./2.*(H(i)+H(i+1));
                end
            end
        end
    end
end
        
        
        
        
        
        
        
        
        