classdef Fig_1
    methods (Static)
        %% > Wrap-up Fig_1.
        function [] = WrapUp_Fig_1(Plot_1,Exp_1,Fig,inp,msh,str,len)
            if Plot_1
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
                %  > Figure.
                figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
                Fig_1.Plot_1(inp,msh,1762,len);
                %  > Export as .pdf.
                if Exp_1
                    Fig_Tools.Export_PDF('Fig_1','../[Figures]/Fig_1');
                end
            end
        end
        
        %% > Auxiliary functions.
        % >> Plot 1.
        function [] = Plot_1(inp,msh,iF,len)
            % >> Local variables.
            Xv_i = inp.msh.lim.Xv_i;
            Xv_f = inp.msh.lim.Xv_f;
            Yv_i = inp.msh.lim.Yv_i;
            Yv_f = inp.msh.lim.Yv_f;
            p    = inp.fr.np;
            NLay = 1./2.*(p+1);
            gq   = false;
            %  > Outer boundary.
            NE   = [Xv_f,Yv_f];
            NW   = [Xv_i,Yv_f];
            SW   = [Xv_i,Yv_i];
            SE   = [Xv_f,Yv_i];
            
            %% > Face selection.
            Sz = nnz(~cellfun(@isempty,msh.s.c(:,iF)));
            C  = linspecer(9,'qualitative');            
            %% > Cell borders.
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            patch([NE(1),NW(1),SW(1),SE(1)],[NE(2),NW(2),SW(2),SE(2)],'w');
            hold on;
            for i = 1:msh.c.NC
                if ismembc(i,msh.f.c{iF})
                    continue;
                end
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),'w');
            end
            % %% > Boundaries.
            % %  > Cells.
            % hold on;
            % for i = 1:size(msh.bnd.c,2)
            %     [f{i}(:,1),f{i}(:,2)] = Fig_Tools.Order_Clockwise('T',msh.bnd.c{1,i});
            %     patch(f{i}(:,1),f{i}(:,2),'r');
            % end
            % %  > Faces.
            % hold on;
            % for i = 1:size(msh.bnd.f,2)
            %     plot(msh.bnd.f{1,i}(:,1),msh.bnd.f{1,i}(:,2),'-b','Linewidth',1.5);
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
            hold on;
            for i = 1:Sz
                if i <= NLay
                    for j = 1:length(msh.s.c{i,iF})
                        [f_1{i,j}(:,1),f_1{i,j}(:,2)] = Fig_Tools.Order_Clockwise('T',msh.c.f.xy_v{msh.s.c{i,iF}(j)});
                        patch(f_1{i,j}(:,1),f_1{i,j}(:,2),C(i,:),'FaceAlpha',0.25,'Linestyle','None');
                    end
                    %  > Cells.
                    c  (i) = plot(msh.c.mean(1,msh.s.c{i,iF}),msh.c.mean(2,msh.s.c{i,iF}),'s','Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',3.5);
                    leg(i) = convertCharsToStrings(num2str(i));
                    %  > Faces.
                    plot(msh.f.mean(1,msh.s.f{i,iF}),msh.f.mean(2,msh.s.f{i,iF}),'s','Color',C(i,:),'MarkerFaceColor',C(i,:),'MarkerSize',3.5);
                else
                    for j = 1:length(msh.s.c{i,iF})
                        [f_1{i,j}(:,1),f_1{i,j}(:,2)] = Fig_Tools.Order_Clockwise('T',msh.c.f.xy_v{msh.s.c{i,iF}(j)});
                        patch(f_1{i,j}(:,1),f_1{i,j}(:,2),C(NLay+1,:),'FaceAlpha',0.25,'Linestyle','None');
                    end
                    %  > Cells.
                    c  (NLay+1) = plot(msh.c.mean(1,msh.s.c{i,iF}),msh.c.mean(2,msh.s.c{i,iF}),'s','Color',C(NLay+1,:),'MarkerFaceColor',C(NLay+1,:),'MarkerSize',3.5);
                    leg(NLay+1) = '$*$';
                    %  > Faces.
                    if ~isempty(msh.s.f{i,iF})
                        plot(msh.f.mean(1,msh.s.f{i,iF}),msh.f.mean(2,msh.s.f{i,iF}),'s','Color',C(NLay+1,:),'MarkerFaceColor',C(NLay+1,:),'MarkerSize',3.5);
                    end
                end
            end
            %  > Face.
            plot(msh.f.xy_v{iF}(:,1),msh.f.xy_v{iF}(:,2),'-','Color',C(1,:),'Linewidth',2.0);
            Fig_Tools.Plot_Limits(inp,msh,iF);
            legend(c,leg,'Interpreter','latex','Location','NortheastOutside','FontSize',12);
            Fig_Tools.ChangeLook_1(inp,len);
        end
    end
end