classdef Fig_3
    methods (Static)         
        %% > Wrap up Fig_3.
        function [] = WrapUp_Fig_3(iD_Var)
            %  > Deal data.
            [n,H,X,H_av,CR] = Fig_3.Process_Data(iD_Var);
            
            % >> Fig_3_1.
            figure(2); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
            for i = 1:3
                %  > Suplots(:,1).
                Fig_3.Subplots_C1(n,i,H,X);
            end
            % >> Fig_3_2.
            figure(3); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
            for i = 1:3
                %  > Suplots(:,2).
                Fig_3.Subplots_C2(n,i,H_av,CR);
            end
        end
         
        %% > Auxiliary functions.
        % >> Process data.
        function [j,H_Sort,X_Sort,H_av,CR] = Process_Data(iD_Var)
            %  > Deal data.
            for i = 1:size(iD_Var,2)
                n(i) = iD_Var{1,i};
            end
            j = unique(n); j = sort(j,'ascend');
            for i = 1:length(j)
                k{i} = find(n == j(i));
                for l = 1:size(k{i},2)
                    H{i}(l) = 1./iD_Var{2,k{i}(l)};
                    X{i}{l} = iD_Var{3,k{i}(l)}.E;
                end
                %  > Sort...
                [~,m{i}] = sort(H{i},'descend');
                for l = 1:size(k{i},2)
                    H_Sort{i}(l) = H{i}(m{i}(l));
                    X_Sort{i}{l} = X{i}{m{i}(l)};
                end
                %  > Compute convergence rate.
                [H_av{i},CR{i}] = Fig_Tools.Compute_ConvergenceRate(H_Sort{i},X_Sort{i});
            end
        end
        
        % >> Subplots.
        %  > Subplots(:,1).
        function [] = Subplots_C1(n,j,H,X)
            %  > Set colors/markers.
            Color = linspecer(size(X,2),'qualitative'); Marker = Fig_Tools.Set_Markers();
            subplot(1,3,j);
            hold on;
            for i = 1:size(X,2)
                for k = 1:size(X{i},2)
                    EX(i,k) = cell2mat(X{i}{k}(j));
                end
                plot(H{i},EX(i,:),Marker{i},'color',Color(i,:),'Linewidth',1.5,'MarkerFaceColor',Color(i,:),'MarkerSize',3.5);
                legendName(i) = convertCharsToStrings(num2str(n(i)));
            end
            for i = 1:size(X,2)
                %  > Get index.
                [~,l] = min(abs(H{i}-1./500));
                fplot(@(x) (EX(i,l)./H{i}(l).^n(i)).*x.^n(i),[H{i}(end),H{i}(1)],'--','color',Color(i,:),'Linewidth',0.50);
                legendName(size(X,2)+i) = convertCharsToStrings(num2str(n(i)));
            end
            legend(legendName,'Interpreter','latex','Location','Northwest','FontSize',10);
            Fig_Tools.ChangeLook_2(j,H);
        end
        %  > Subplots(:,2).
        function [] = Subplots_C2(n,j,H,X)
            %  > Set colors/markers.
            Color = linspecer(size(X,2),'qualitative'); Marker = Fig_Tools.Set_Markers();
            subplot(1,3,j);
            hold on;
            for i = 1:size(X,2)
                plot(H{i},X{i}(:,j),Marker{i},'color',Color(i,:),'Linewidth',1.5,'MarkerFaceColor',Color(i,:),'MarkerSize',3.5);
                legendName(i) = convertCharsToStrings(num2str(n(i)));
            end
            legend(legendName,'Interpreter','latex','Location','Northwest','FontSize',10);
            Fig_Tools.ChangeLook_2(j,H);
        end
    end
end    
      
        
        
        