classdef Fig_1
    methods (Static)         
        %% > Wrap up Fig_1.
        function [] = WrapUp_Fig_1(Fig,inp,msh)
            %  > iD: Select random bulk/bnd cell.
            str = 'blk';
            bnd = cell2mat(msh.bnd.c(2,:));
            if strcmpi(str,'bnd')
                iD = bnd(randperm(length(bnd),1));
            elseif strcmpi(str,'blk')
                blk = setdiff(1:msh.c.NC,bnd);
                if isempty(blk)
                    return;
                end
                iD  = blk(randperm(length(blk),1));
            end
            
            %  > Figure.
            figure(Fig(1)); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
            Fig_1.Plot_1(inp,msh,iD,0);
        end
        
        %% > Auxiliary functions.
        % >> Plot 1.
        function [] = Plot_1(inp,msh,iD,len)
            %% > Face selection.
            C  = linspecer(size(msh.s.st,1),'qualitative');
            %  > Identify face(iD,fD).
            Flag = zeros(1,size(msh.f.nb,2));
            for i = 1:size(msh.f.nb,2)
                if ismembc(iD,msh.f.nb{i})
                    Flag(i) = i;
                end
            end
            fi = find(Flag);
            fD = fi(randperm(length(fi),1));
        
            %% > Cell borders.
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                patch(msh.c.XY_v{i}(:,1),msh.c.XY_v{i}(:,2),'w');
            end
%             %% > Boundaries.
%             %  > Cells.
%             hold on;
%             for i = 1:size(msh.bnd.c,2)
%                 %  > Test 1.
%                 f{i} = unique(cell2mat(reshape(msh.bnd.c{1,i},[size(msh.bnd.c{1,i},2),1])),'rows','stable'); [f{i}(:,1),f{i}(:,2)] = Fig_Tools.Order_Clockwise(f{i}(:,1),f{i}(:,2));
%                 patch(f{i}(:,1),f{i}(:,2),'r','FaceAlpha',0.25);
%                 %  > Test 2.
%                 %  patch(msh.c.XY_v{msh.bnd.c{2,i}}(:,1),msh.c.XY_v{msh.bnd.c{2,i}}(:,2),'b','FaceAlpha',0.25);
%             end
%             %  > Faces.
%             hold on;
%             for i = 1:size(msh.bnd.f,2)
%                 %  > Test 1.
%                 %  plot(msh.bnd.f{1,i}(:,1),msh.bnd.f{1,i}(:,2),'-r','Linewidth',2.5);
%                 %  > Test 2.
%                 plot(msh.f.XY_v{msh.bnd.f{2,i}}(:,1),msh.f.XY_v{msh.bnd.f{2,i}}(:,2),'-r','Linewidth',2.5);
%             end
            %% > Centroids.
            %  > Cells.
            hold on;
            for i = 1:msh.c.NC
                plot(msh.c.mean(1,i),msh.c.mean(2,i),'ok','MarkerFaceColor','k','MarkerSize',2.5);
            end
%             %  > Faces.
%             hold on;
%             for j = 1:msh.f.NF
%                 plot(msh.f.mean(1,j),msh.f.mean(2,j),'ob','MarkerFaceColor','b','MarkerSize',2.5);
%             end
%             %% > Normals.
%             hold on;
%             for i = 1:msh.c.NC
%                 for j = 1:size(msh.c.f.XY_v{i},2)
%                     quiver(msh.c.f.mean{i}(1,j),msh.c.f.mean{i}(2,j),msh.c.f.Nf{i}(1,j).*len,msh.c.f.Nf{i}(2,j).*len,'AutoScale','off','Color','b');
%                 end
%             end
%             %% > Neighbours(iD).
%             %  > Cell(iD).
%             hold on;
%             for i = 1:length(msh.c.nb{iD})
%                 [f{i}(:,1),f{i}(:,2)] = Fig_Tools.Order_Clockwise(msh.c.XY_v{msh.c.nb{iD}(i)}(:,1),msh.c.XY_v{msh.c.nb{iD}(i)}(:,2));
%                 patch(f{i}(:,1),f{i}(:,2),'b','FaceAlpha',0.25);
%             end
%             %  > Face(iD).
%             hold on;
%             plot(msh.f.XY_v{iD}(:,1),msh.f.XY_v{iD}(:,2),'-b','Linewidth',2.5);
%             for i = 1:length(msh.f.nb{iD})
%                 f{i} = unique(cell2mat(reshape(msh.c.f.XY_v{msh.f.nb{iD}(i)},[size(msh.c.f.XY_v{msh.f.nb{iD}(i)},2),1])),'rows','stable'); [f{i}(:,1),f{i}(:,2)] = Fig_Tools.Order_Clockwise(f{i}(:,1),f{i}(:,2));
%                 patch(f{i}(:,1),f{i}(:,2),'b','FaceAlpha',0.25);
%             end
            %% > Stencil.
            % >> Stencil representation.
            hold on;
            plot(msh.f.XY_v{fD}(:,1),msh.f.XY_v{fD}(:,2),'-b','Linewidth',2.0);
            for i = 1:size(msh.s.st,1)
                %  > Cells.
                for j = 1:length(msh.s.st{i,fD})
                    %  > To patch...
                    f{i,j} = unique(cell2mat(reshape(msh.c.f.XY_v{msh.s.st{i,fD}(j)},[size(msh.c.f.XY_v{msh.s.st{i,fD}(j)},2),1])),'rows','stable'); [f{i,j}(:,1),f{i,j}(:,2)] = Fig_Tools.Order_Clockwise(f{i,j}(:,1),f{i,j}(:,2));
                    %  > Pacth cell(s).
                    patch(f{i,j}(:,1),f{i,j}(:,2),C(i,:),'FaceAlpha',0.25);
                end
                %  > Points.
                for j = 1:length(msh.s.XY_v{i,fD})
                    p_leg(i) = plot(msh.s.XY_v{i,fD}(1,j),msh.s.XY_v{i,fD}(2,j),'s','Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',5.0);
                end
                leg(i) = convertCharsToStrings(num2str(i));
                %  > Gauss integration points.
                for j = 1:size(msh.f.fg{fD}.Points,1)
                    scatter(msh.f.fg{fD}.Points(j,1),msh.f.fg{fD}.Points(j,2),50.*msh.f.fg{fD}.Weights(j),'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r');
                end
            end
            legend(p_leg,leg,'Interpreter','latex','Location','NortheastOutside','FontSize',10);
            Fig_Tools.ChangeLook_1(inp,len); 
        end
    end
end