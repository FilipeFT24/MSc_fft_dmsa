classdef Test_1
    methods (Static)
        %% > Wrap-up Test_1.
        function [] = WrapUp_T1(Plot_T1,Exp_T1,Fig,X)
            if Plot_T1
                %  > Figure #1.
                figure(Fig(1)); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
                Test_1.Plot_1(X,1);
                Test_1.Plot_1(X,2);
                Test_1.Plot_1(X,3);
                figure(Fig(2)); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
                Test_1.Plot_2(X,1);
                Test_1.Plot_2(X,2);
                Test_1.Plot_2(X,3);
                %  > Figure #2.
                %  > Export as .pdf.
                if Exp_T1
                    Fig_Tools.Export_PDF('Test_1','../[Figures]/Test_1');
                end
            end
        end
        
        %% > Auxiliary functions.
        % >> Plot 1.
        function [] = Plot_1(X,j)
            %  > Set colors/markers.
            Color = linspecer(size(X.NR.E,2),'qualitative'); Marker = Fig_Tools.Set_Markers();
            subplot(1,3,j);
            hold on;
            for i = 1:size(X.NR.E,2)
                for k = 1:size(X.NR.E{i},2)
                    Ej(i,k) = X.NR.E{i}{k}(j);
                end
                plot(X.NR.H{i},Ej(i,:),Marker{i},'color',Color(i,:),'Linewidth',1.5,'MarkerFaceColor',Color(i,:),'MarkerSize',3.5);
                legendName(i) = convertCharsToStrings(num2str(X.n(i)));
            end
            for i = 1:size(X.NR.E,2)
                %  > Get index...
                [~,l] = min(abs(X.NR.H{i}-5.0E-02));
                fplot(@(x) (Ej(i,l)./X.NR.H{i}(l).^X.n(i)).*x.^X.n(i),[X.NR.H{i}(end),X.NR.H{i}(1)],'--','color',Color(i,:),'Linewidth',0.50);
                legendName(size(X.NR.E,2)+i) = convertCharsToStrings(num2str(X.n(i)));
            end
            legend(legendName,'Interpreter','latex','Location','Northeast','FontSize',10);
            Fig_Tools.ChangeLook_2(j,X.NR.H);
        end
        % >> Plot 2.
        function [] = Plot_2(X,j) 
            %  > Set colors/markers.
            Color = linspecer(size(X.NR.E,2),'qualitative'); Marker = Fig_Tools.Set_Markers();
            subplot(1,3,j);
            hold on;
            for i = 1:size(X.NR.E,2)
                plot(X.CR.H{i},X.CR.R{i}(:,j),Marker{i},'color',Color(i,:),'Linewidth',1.5,'MarkerFaceColor',Color(i,:),'MarkerSize',3.5);
                legendName(i) = convertCharsToStrings(num2str(X.n(i)));
            end
            legend(legendName,'Interpreter','latex','Location','Northwest','FontSize',10);
            Fig_Tools.ChangeLook_2(j,X.NR.H);
        end
    end
end