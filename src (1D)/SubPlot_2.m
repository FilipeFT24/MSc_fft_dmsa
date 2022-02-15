
            %  > Method's order.
            for j = 1:size(msh.s.xf,1)
                for i = 1:msh.f.NF
                    l(j,i) = length(msh.s.xf{j,i});
                end
                k_min(j) = min(l(j,:));
                k_max(j) = max(l(j,:));
            end
%  > 1.1.2. ---------------------------------------------------
            %subplot(7,1,6);
            %Fig_2_1D.SubPlot_2(C,min(k_min),max(k_max));
            %  > 1.1.3. ---------------------------------------------------
            %subplot(7,1,7);
            %Fig_2_1D.SubPlot_3(msh,C,l);
%  > 1.1.2. -------------------------------------------------------
        function [] = SubPlot_2(C,k_min,k_max)
            %  > Legend.
            j = 1;
            hold on;
            for i = k_min:k_max
                if rem(i,2) == 1
                    PN(j) = plot(NaN,NaN,'o','Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',4.5);
                else
                    PN(j) = plot(NaN,NaN,'s','Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',5.0);
                end
                LN{j} = num2str(i);
                j     = j+1;
            end
            legend(PN,LN,'Interpreter','latex','Location','Northeast','FontSize',10,'NumColumns',1);
            %  > Axis.
            ax = gca; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
        end
        %  > 1.1.3. -------------------------------------------------------
        function [] = SubPlot_3(msh,C,l)
            hold on;
            plot(msh.f.Xv,repelem(0.00,msh.f.NF),'-','Color','k','MarkerFaceColor','k','MarkerSize',10.0,'Linewidth',0.5);
            plot(msh.f.Xv,repelem(0.05,msh.f.NF),'-','Color','k','MarkerFaceColor','k','MarkerSize',10.0,'Linewidth',0.5);
            y = [0.05,0];
            for j = 1:size(l,1)
                for i = 1:msh.f.NF
                    if rem(l(j,i),2) == 1
                        plot(msh.f.Xv(i),y(j),'o','Color',C(l(j,i),:),'MarkerFaceColor',C(l(j,i),:),'MarkerSize',4.5);
                    else
                        plot(msh.f.Xv(i),y(j),'s','Color',C(l(j,i),:),'MarkerFaceColor',C(l(j,i),:),'MarkerSize',5.0);
                    end
                end
            end
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(false,true,false,msh.f.Xv,10,"$x$","$y$",20,12);
            ax = gca; ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
        end