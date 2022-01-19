classdef Fig_3
    methods (Static)
        %% > Wrap up Fig_3.
        function [] = WrapUp_Fig_3(Fig,inp,msh,pde,len)
            %  > Figure.
            figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
            %  > #1.
            subplot(1,3,1);
            Fig_3.Plot_1(inp,msh,pde,len);
            %  > #2.
            subplot(1,3,2);
            Fig_3.Plot_2(inp,msh,pde,len);
            %  > #3.
            subplot(1,3,3);
            Fig_3.Plot_3(inp,msh,pde,len);
        end
        
        %% > Auxiliary functions.
        % >> Plot 1.
        function [] = Plot_1(inp,msh,pde,len)
            %% > Cell borders.
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                patch(msh.c.xy_v{1,i}(:,1),msh.c.xy_v{1,i}(:,2),pde.blk.f(i));
            end
            c = Fig_Tools.Colormap_cmocean('thermal');
            c.Label.String = '$\phi$';
            Fig_Tools.ChangeLook_1(inp,len);
            AdvancedColormap('thermal');
            % >> Options:
            %  > AdvancedColormap('cool');
            %  > AdvancedColormap('hot');
            %  > AdvancedColormap('hsv');
            %  > AdvancedColormap('spring');
        end
        % >> Plot 2.
        function [] = Plot_2(inp,msh,pde,len)
            %% > Cell borders.
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                patch(msh.c.xy_v{1,i}(:,1),msh.c.xy_v{1,i}(:,2),pde.Phi(i));
            end
            c = Fig_Tools.Colormap_cmocean('thermal');
            c.Label.String = '$\phi_{_{\mathrm{PDE}}}$';
            Fig_Tools.ChangeLook_1(inp,len);
            AdvancedColormap('thermal');
        end  
        % >> Plot 3.
        function [] = Plot_3(inp,msh,pde,len)
            %% > Cell borders.
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                patch(msh.c.xy_v{1,i}(:,1),msh.c.xy_v{1,i}(:,2),abs(pde.blk.f(i)-pde.Phi(i)));
            end
            % for i = 1:msh.c.NC
            %     scatter(pde.c.Qp{i}.Points(:,1),pde.c.Qp{i}.Points(:,2),10.*pde.c.Qp{i}.Weights./sum(pde.c.Qp{i}.Weights),'Marker','o','MarkerEdgeColor','r','MarkerFaceColor','r');
            % end
            c = Fig_Tools.Colormap_cmocean('thermal');
            c.Label.String = '$\textrm{Absolute error}\left(\epsilon_{abs}\right)$';
            Fig_Tools.ChangeLook_1(inp,len);
            AdvancedColormap('thermal');
        end
    end
end