classdef Fig_4
    methods (Static)
        %% > Wrap up-Fig_4.
        function [] = WrapUp_Fig_4(Plot_4,Fig,inp,msh,pde,len)
            if Plot_4
                %  > Figure.
                figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
                Fig_4.Plot_1(inp,msh,pde,len);
            end
        end
        
        %% > Auxiliary functions.
        % >> Plot 1.
        function [] = Plot_1(inp,msh,pde,len)
            %% > Cell borders.
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                patch(msh.c.xy_v{1,i}(:,1),msh.c.xy_v{1,i}(:,2),abs(pde.blk.f(i)-pde.Phi(i)));
            end
            c = Fig_Tools.Colormap_cmocean('thermal');
            c.Label.String = '$\textrm{Absolute error}\left(\epsilon_{abs}\right)$';
            Fig_Tools.ChangeLook_1(inp,len);
            AdvancedColormap('thermal');
        end
    end
end