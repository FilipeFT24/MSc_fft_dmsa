classdef Fig_1
    methods (Static)
        %% > Wrap up Fig_1.
        function [] = WrapUp_Fig_1(Fig,inp,msh,str,len)
            %  > Auxiliary arrays.
            bnd_f = cell2mat(msh.bnd.f(2,:));
            bnd_c = cell2mat(msh.bnd.f(3,:));
            %  > Select...
            if strcmpi(str,'bnd')
                iX  = randperm(size(msh.bnd.f,2),1);
                iF  = bnd_f(iX);
            elseif strcmpi(str,'blk')
                blk = setdiff(1:msh.f.NF,bnd_f);
                iX  = randperm(length(blk),1);
                iF  = blk(iX);
            end
            figure(Fig(1)); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
            Fig_1.Plot_1(inp,msh,1,len);
        end
        
        %% > Auxiliary functions.
        % >> Plot 1.
        function [] = Plot_1(inp,msh,iF,len)
            %% > Face selection.
            Sz = nnz(~cellfun(@isempty,msh.s.c(:,iF)));
            C  = linspecer(9,'qualitative');
            
            %% > Cell borders.
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),'w');
            end
            % %% > Boundaries.
            % %  > Cells.
            % hold on;
            % for i = 1:size(msh.bnd.c,2)
            %     [f{i}(:,1),f{i}(:,2)] = Fig_Tools.Order_Clockwise('T',msh.bnd.c{1,i});
            %     patch(f{i}(:,1),f{i}(:,2),'r','FaceAlpha',0.10);
            % end
            % %  > Faces.
            % hold on;
            % for i = 1:size(msh.bnd.f,2)
            %     plot(msh.bnd.f{1,i}(:,1),msh.bnd.f{1,i}(:,2),'-r','Linewidth',2.5);
            % end
            %% > Centroids.
            %  > Cells.
            hold on;
            for i = 1:msh.c.NC
                plot(msh.c.mean(1,i),msh.c.mean(2,i),'ok','MarkerFaceColor','k','MarkerSize',1.5);
            end
            % %  > Faces.
            % hold on;
            % for j = 1:msh.f.NF
            %     plot(msh.f.mean(1,j),msh.f.mean(2,j),'ob','MarkerFaceColor','b','MarkerSize',1.5);
            % end
            % %% > Normals.
            % hold on;
            % for i = 1:msh.c.NC
            %     for j = 1:size(msh.c.f.xy_v{i},2)
            %         quiver(msh.c.f.mean{i}(1,j),msh.c.f.mean{i}(2,j),msh.c.f.Nf{i}(1,j).*len,msh.c.f.Nf{i}(2,j).*len,'AutoScale','off','Color','b');
            %     end
            % end
            % %% > Neighbours(iF).
            % %  > Cells.
            % hold on;
            % for i = 1:length(msh.f.c{iF})
            %     [f{i}(:,1),f{i}(:,2)] = Fig_Tools.Order_Clockwise('F',msh.c.xy_v{msh.f.c{iF}(i)});
            %     patch(f{i}(:,1),f{i}(:,2),'b','FaceAlpha',0.10);
            % end
            %% > Stencil.
            p    = inp.fr.np;
            NLay = 1./2.*(p+1);
            hold on;
            for i = 1:Sz
                if i <= NLay
                    for j = 1:length(msh.s.c{i,iF})
                        [f_1{i,j}(:,1),f_1{i,j}(:,2)] = Fig_Tools.Order_Clockwise('T',msh.c.f.xy_v{msh.s.c{i,iF}(j)});
                        patch(f_1{i,j}(:,1),f_1{i,j}(:,2),C(i,:),'FaceAlpha',0.25);
                    end
                    %  > Cells.
                    c  (i) = plot(msh.c.mean(1,msh.s.c{i,iF}),msh.c.mean(2,msh.s.c{i,iF}),'s','Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',3.5);
                    leg(i) = convertCharsToStrings(num2str(i));
                    %  > Faces.
                    plot(msh.f.mean(1,msh.s.f{i,iF}),msh.f.mean(2,msh.s.f{i,iF}),'s','Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',3.5);                
                else
                    for j = 1:length(msh.s.c{i,iF})
                        [f_1{i,j}(:,1),f_1{i,j}(:,2)] = Fig_Tools.Order_Clockwise('T',msh.c.f.xy_v{msh.s.c{i,iF}(j)});
                        patch(f_1{i,j}(:,1),f_1{i,j}(:,2),C(NLay+1,:),'FaceAlpha',0.25);
                    end
                    %  > Cells.
                    c  (NLay+1) = plot(msh.c.mean(1,msh.s.c{i,iF}),msh.c.mean(2,msh.s.c{i,iF}),'s','Color',C(NLay+1,:),'MarkerFaceColor',C(NLay+1,:),'MarkerSize',3.5);
                    leg(NLay+1) = '*';
                    %  > Faces.
                    %if ~isempty(msh.s.f{i,iF})
                    %    plot(msh.f.mean(1,msh.s.f{i,iF}),msh.f.mean(2,msh.s.f{i,iF}),'s','Color',C(NLay+1,:),'MarkerFaceColor',C(NLay+1,:),'MarkerSize',3.5);
                    %end
                end
            end
            
            
            for i = 6:6
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),'m');
            end
            for i = 14:14
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),'m');
            end
            for i = 22:22
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),'m');
            end
            for i = 30:30
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),'m');
            end
            for i = 38:38
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),'m');
            end
%             for i = 41:45
%                 patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),'m');
%             end
            
            
            
            %  > Face.
            plot(msh.f.xy_v{iF}(:,1),msh.f.xy_v{iF}(:,2),'-','Color',C(1,:),'Linewidth',2.5);
            %  > Limits.
            Fig_Tools.Plot_Limits(inp,msh,iF);
            %  > Other stuff.
            legend(c,leg,'Interpreter','latex','Location','NortheastOutside','FontSize',10);
            Fig_Tools.ChangeLook_1(inp,len);
        end
    end
end