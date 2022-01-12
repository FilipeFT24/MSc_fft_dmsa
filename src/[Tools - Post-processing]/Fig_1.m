classdef Fig_1
    methods (Static)         
        %% > Wrap up Fig_1.
        function [] = WrapUp_Fig_1(Fig,inp,msh)
            %  > iD: Select random bulk/bnd cell.
            str = 'bnd';
            bnd = cell2mat(msh.bnd.c(2,:));
            if strcmpi(str,'bnd')
                iC = bnd(randperm(length(bnd),1));
            elseif strcmpi(str,'blk')
                blk = setdiff(1:msh.c.NC,bnd);
                if isempty(blk)
                    return;
                end
                iC  = blk(randperm(length(blk),1));
            end
            
            %  > Figure.
            figure(Fig(1)); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
            Fig_1.Plot_1(inp,msh,iC,0);
        end
        
        %% > Auxiliary functions.
        % >> Plot 1.
        function [] = Plot_1(inp,msh,iD,len)
            %% > Face selection.
            C  = linspecer(size(msh.s.st,1),'qualitative');
            fD = msh.c.f.faces{iD}(randperm(length(msh.c.f.faces{iD}),1));
        
            %% > Cell borders.
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),'w');
            end
%             %% > Boundaries.
%             %  > Cells.
%             hold on;
%             for i = 1:size(msh.bnd.c,2)
%                 %  > Test 1.
%                 f{i} = unique(cell2mat(reshape(msh.bnd.c{1,i},[size(msh.bnd.c{1,i},2),1])),'rows','stable'); [f{i}(:,1),f{i}(:,2)] = Fig_Tools.Order_Clockwise(f{i}(:,1),f{i}(:,2));
%                 patch(f{i}(:,1),f{i}(:,2),'r','FaceAlpha',0.25);
%                 %  > Test 2.
%                 %  patch(msh.c.xy_v{msh.bnd.c{2,i}}(:,1),msh.c.xy_v{msh.bnd.c{2,i}}(:,2),'b','FaceAlpha',0.25);
%             end
%             %  > Faces.
%             hold on;
%             for i = 1:size(msh.bnd.f,2)
%                 %  > Test 1.
%                 %  plot(msh.bnd.f{1,i}(:,1),msh.bnd.f{1,i}(:,2),'-r','Linewidth',2.5);
%                 %  > Test 2.
%                 plot(msh.f.xy_v{msh.bnd.f{2,i}}(:,1),msh.f.xy_v{msh.bnd.f{2,i}}(:,2),'-r','Linewidth',2.5);
%             end
            %% > Centroids.
            %  > Cells.
            hold on;
            for i = 1:msh.c.NC
                plot(msh.c.mean(1,i),msh.c.mean(2,i),'ok','MarkerFaceColor','k','MarkerSize',1.5);
            end
%             %  > Faces.
%             hold on;
%             for j = 1:msh.f.NF
%                 plot(msh.f.mean(1,j),msh.f.mean(2,j),'ob','MarkerFaceColor','b','MarkerSize',2.5);
%             end
%             %% > Normals.
%             hold on;
%             for i = 1:msh.c.NC
%                 for j = 1:size(msh.c.f.xy_v{i},2)
%                     quiver(msh.c.f.mean{i}(1,j),msh.c.f.mean{i}(2,j),msh.c.f.Nf{i}(1,j).*len,msh.c.f.Nf{i}(2,j).*len,'AutoScale','off','Color','b');
%                 end
%             end
%             %% > Neighbours(iD).
%             %  > Cell(iD).
%             hold on;
%             for i = 1:length(msh.c.nb{iD})
%                 [f{i}(:,1),f{i}(:,2)] = Fig_Tools.Order_Clockwise(msh.c.xy_v{msh.c.nb{iD}(i)}(:,1),msh.c.xy_v{msh.c.nb{iD}(i)}(:,2));
%                 patch(f{i}(:,1),f{i}(:,2),'b','FaceAlpha',0.25);
%             end
            %% > Stencil.
            % >> Stencil representation.
            hold on;
            plot(msh.f.xy_v{fD}(:,1),msh.f.xy_v{fD}(:,2),'-b','Linewidth',2.0);
            for i = 1:size(msh.s.st,1)
                %  > Cells.
                for j = 1:length(msh.s.st{i,fD})
                    %  > To patch...
                    f{i,j} = unique(cell2mat(reshape(msh.c.f.xy_v{msh.s.st{i,fD}(j)},[size(msh.c.f.xy_v{msh.s.st{i,fD}(j)},2),1])),'rows','stable'); [f{i,j}(:,1),f{i,j}(:,2)] = Fig_Tools.Order_Clockwise(f{i,j}(:,1),f{i,j}(:,2));
                    %  > Pacth cell(s).
                    patch(f{i,j}(:,1),f{i,j}(:,2),C(i,:),'FaceAlpha',0.25);
                end
                %  > Points.
                for j = 1:length(msh.s.xy_v{i,fD})
                    p_leg(i) = plot(msh.s.xy_v{i,fD}(1,j),msh.s.xy_v{i,fD}(2,j),'s','Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',3.5);
                end
                leg(i) = convertCharsToStrings(num2str(i));
                %  > Gauss integration points.
                for j = 1:size(msh.f.gq{fD}.Points,1)
                    scatter(msh.f.gq{fD}.Points(j,1),msh.f.gq{fD}.Points(j,2),50.*msh.f.gq{fD}.Weights(j),'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r');
                end
            end
            legend(p_leg,leg,'Interpreter','latex','Location','NortheastOutside','FontSize',10);
            Fig_Tools.ChangeLook_1(inp,len); 
        end
    end
end