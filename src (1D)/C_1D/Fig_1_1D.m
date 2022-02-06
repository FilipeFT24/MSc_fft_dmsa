classdef Fig_1_1D
    methods (Static)         
        %% > Wrap-up Fig_1 (1D).
        function [] = WrapUp_Fig_1_1D(Plot_1,Exp_1,Fig,msh)
            if Plot_1
                %  > Select...
                non_empty = find(cellfun(@isempty,msh.s.f));
                iF        = non_empty(randperm(length(non_empty),1));               
                %  > Figure.
                figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
                Fig_1_1D.Plot(msh,iF);
                %  > Export as .pdf.
                if Exp_1
                    Fig_Tools_1D.Export_PDF('Fig_1','../[Figures]/[1D]/Fig_1');
                end
            end
        end
        
        %% > Auxiliary functions.
        function [] = Plot(msh,iF)
            %  > Patch...
            wdt  = 0.05;
            c_xy = Fig_Tools_1D.ToPatch(msh,wdt);
            C    = linspecer(9,'qualitative');            
            %% > #1.
            subplot(3,1,1);
            %  > Organize layers...
            l1    = length(msh.s.x_v_t{1});
            L1    = msh.s.c{1};
            L1    = [0,L1];
            lay_1 = reshape(L1,2,l1./2)';
            %  > Plot...
            hold on;
            for i = 1:msh.c.NC
                patch(c_xy{i}(1,:),c_xy{i}(2,:),'w');
            end
            hold on;
            for i = 1:l1./2
                if i == 1
                    for j = 2
                        patch(c_xy{lay_1(i,j)}(1,:),c_xy{lay_1(i,j)}(2,:),C(i,:),'FaceAlpha',0.25,'Linestyle','None');
                    end
                    c_1(i) = plot([msh.f.Xv(1),msh.c.Xc(lay_1(1,size(lay_1,2)))],[0,0],'s','Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',4.75);
                else
                    for j = 1:2
                        patch(c_xy{lay_1(i,j)}(1,:),c_xy{lay_1(i,j)}(2,:),C(i,:),'FaceAlpha',0.25,'Linestyle','None');
                    end
                    c_1(i) = plot(msh.c.Xc(lay_1(i,:)),[0,0],'s','Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',4.75);
                end
                leg_1(i) = convertCharsToStrings(num2str(i));
            end
            plot([msh.f.Xv(1),msh.f.Xv(1)],[-wdt,wdt],'-','Color',C(1,:),'Linewidth',2.0);
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(false,msh.f.Xv,10,30);
            ax = gca; 
            set(gca,'XTick',[]); ax.YAxis.Visible = 'off';
            legend(c_1,leg_1,'Interpreter','latex','Location','NortheastOutside','FontSize',12);     
            %% > #2.
            subplot(3,1,2);
            %  > Organize layers...
            lF   = length(msh.s.x_v_t{msh.f.NF});
            LF   = msh.s.c{msh.f.NF};
            LF   = [LF,0];
            layF = flipud(reshape(LF,2,lF./2)');
            %  > Plot...
            hold on;
            for i = 1:msh.c.NC
                patch(c_xy{i}(1,:),c_xy{i}(2,:),'w');
            end
            hold on;
            for i = 1:lF./2
                if i == 1
                    for j = 1
                        patch(c_xy{layF(i,j)}(1,:),c_xy{layF(i,j)}(2,:),C(i,:),'FaceAlpha',0.25,'Linestyle','None');
                    end
                    c_F(i) = plot([msh.f.Xv(msh.f.NF),msh.c.Xc(layF(1,1))],[0,0],'s','Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',4.75);
                else
                    for j = 1:2
                        patch(c_xy{layF(i,j)}(1,:),c_xy{layF(i,j)}(2,:),C(i,:),'FaceAlpha',0.25,'Linestyle','None');
                    end
                    c_F(i) = plot(msh.c.Xc(layF(i,:)),[0,0],'s','Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',4.75);
                end
                leg_F(i) = convertCharsToStrings(num2str(i));
            end
            plot([msh.f.Xv(msh.f.NF),msh.f.Xv(msh.f.NF)],[-wdt,wdt],'-','Color',C(1,:),'Linewidth',2.0);
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(false,msh.f.Xv,10,30);
            ax = gca; 
            set(gca,'XTick',[]); ax.YAxis.Visible = 'off';
            legend(c_F,leg_F,'Interpreter','latex','Location','NortheastOutside','FontSize',12);         
            %% > #3.
            subplot(3,1,3);
            %  > Organize layers...
            li = length(msh.s.x_v_t{iF});
            Li = msh.s.c{iF};
            j  = 1:2;
            for i = 1:li./2
                lay_i(i,j) = [Li(li./2+1-i),Li(li./2+i)];
            end
            %  > Plot...
            hold on;
            for i = 1:msh.c.NC
                patch(c_xy{i}(1,:),c_xy{i}(2,:),'w');
            end
            hold on;
            for i = 1:li./2
                for k = j
                    patch(c_xy{lay_i(i,k)}(1,:),c_xy{lay_i(i,k)}(2,:),C(i,:),'FaceAlpha',0.25,'Linestyle','None');
                end
                c_i  (i) = plot(msh.c.Xc(lay_i(i,:)),[0,0],'s','Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',4.75);
                leg_i(i) = convertCharsToStrings(num2str(i));
            end
            plot([msh.f.Xv(iF),msh.f.Xv(iF)],[-wdt,wdt],'-','Color',C(1,:),'Linewidth',2.0);
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(false,msh.f.Xv,10,30);
            ax = gca; 
            set(gca,'XTick',[]); ax.YAxis.Visible = 'off';
            legend(c_i,leg_i,'Interpreter','latex','Location','NortheastOutside','FontSize',12);            
        end
    end
end      