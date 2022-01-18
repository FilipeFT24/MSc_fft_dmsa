classdef Fig_2
    methods (Static)
        %% > Wrap up Fig_2.
        function [] = WrapUp_Fig_2(Fig,inp,msh,pde,len)
            figure(Fig(1)); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
            %  > #1.
            subplot(1,2,1);
            Fig_2.Plot_1(inp,msh,pde,len);
            %  > #2.
            subplot(1,2,2);
            Fig_2.Plot_2(inp,msh,pde,len);
        end
        
        %% > Auxiliary functions.
        % >> Plot 1.
        function [] = Plot_1(inp,msh,pde,len)           
            %% > Condition number (Df).
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                c_Df(i) = mean(pde.mat.cd_Df(msh.c.f.f{i}));
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),c_Df(i));
            end
            %  > Mean condition number.
            mean_c_Df = mean(pde.mat.cd_Df);
            %  > Index of maximum condition number.
            for i = 1:msh.f.NF
                if pde.mat.cd_Df(i) > mean_c_Df
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
                c_Dwf(i) = mean(pde.mat.cd_Dwf(msh.c.f.f{i}));
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),c_Dwf(i));
            end
            %  > Mean condition number.
            mean_c_Dwf = mean(pde.mat.cd_Dwf);
            %  > Index of maximum condition number.
            for i = 1:msh.f.NF
                if pde.mat.cd_Dwf(i) > mean_c_Dwf
                    plot(msh.f.xy_v{i}(:,1),msh.f.xy_v{i}(:,2),'-r','Linewidth',2.0);
                end
            end
            
            %  > Other stuff.
            c = Fig_Tools.Colormap_cmocean('thermal');
            c.Label.String = '$\textrm{cond}\left(D_{wf}\right)$';
            Fig_Tools.ChangeLook_1(inp,len);
            AdvancedColormap('thermal');           
        end
    end
end