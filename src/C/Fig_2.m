classdef Fig_2
    methods (Static)
        %% > Wrap up Fig_2.
        function [] = WrapUp_Fig_2(Fig,inp,msh,pde,len)
            %  > Figure.
            figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
            %  > #1.
            subplot(1,3,1);
            Fig_2.Plot_1(inp,msh,pde,len);
            %  > #2.
            subplot(1,3,2);
            Fig_2.Plot_2(inp,msh,pde,len);
            %  > #3.
            subplot(1,3,3);
            Fig_2.Plot_3(inp,msh,pde,len);
        end
        
        %% > Auxiliary functions.
        % >> Plot 1.
        function [] = Plot_1(inp,msh,pde,len)           
            %% > Condition number (Df).
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                for j = 1:length(msh.c.f.f{i})
                    c_Df_ij(i,j) = cond(pde.f.Df{msh.c.f.f{i}(j)});
                end
                c_Df_i(i) = mean(c_Df_ij(i,:));
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),c_Df_i(i));
            end
            %  > Mean condition number.
            mean_c_Df_i = mean(c_Df_i);
            %  > Index of maximum condition number.
            for i = 1:msh.f.NF
                if cond(pde.f.Df{i}) > mean_c_Df_i
                    plot(msh.f.xy_v{i}(:,1),msh.f.xy_v{i}(:,2),'-r','Linewidth',2.0);
                end
            end

            %  > Other stuff.
            c = Fig_Tools.Colormap_cmocean('thermal');
            c.Label.String = '$\textrm{cond}\left(D_{f}\right)$';
            Fig_Tools.ChangeLook_1(inp,len);
            AdvancedColormap('thermal');           
        end
        % >> Plot 2.
        function [] = Plot_2(inp,msh,pde,len)           
            %% > Condition number (Dwf).
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                for j = 1:length(msh.c.f.f{i})
                    c_Dfw_ij(i,j) = cond(pde.f.Dwf{msh.c.f.f{i}(j)});
                end
                c_Dfw_i(i) = mean(c_Dfw_ij(i,:));
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),c_Dfw_i(i));
            end
            %  > Mean condition number.
            mean_c_Dwf_i = mean(c_Dfw_i);
            %  > Index of maximum condition number.
            for i = 1:msh.f.NF
                if cond(pde.f.Dwf{i}) > mean_c_Dwf_i
                    plot(msh.f.xy_v{i}(:,1),msh.f.xy_v{i}(:,2),'-r','Linewidth',2.0);
                end
            end
            c = Fig_Tools.Colormap_cmocean('thermal');
            c.Label.String = '$\textrm{cond}\left(D_{wf}\right)$';
            Fig_Tools.ChangeLook_1(inp,len);
            AdvancedColormap('thermal');           
        end
        % >> Plot 3.
        function [] = Plot_3(inp,msh,pde,len)           
            %% > Condition number (Pf).
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                for j = 1:length(msh.c.f.f{i})
                    c_Pf_ij(i,j) = cond(pde.f.Pf{msh.c.f.f{i}(j)});
                end
                c_Pf_i(i) = mean(c_Pf_ij(i,:));
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),c_Pf_i(i));
            end
            %  > Mean condition number.
            mean_c_Pf_i = mean(c_Pf_i);
            %  > Index of maximum condition number.
            for i = 1:msh.f.NF
                if cond(pde.f.Pf{i}) > mean_c_Pf_i
                    plot(msh.f.xy_v{i}(:,1),msh.f.xy_v{i}(:,2),'-r','Linewidth',2.0);
                end
            end
            c = Fig_Tools.Colormap_cmocean('thermal');
            c.Label.String = '$\textrm{cond}\left(P_{f}\right)$';
            Fig_Tools.ChangeLook_1(inp,len);
            AdvancedColormap('thermal');           
        end
    end
end