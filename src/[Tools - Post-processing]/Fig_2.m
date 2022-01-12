classdef Fig_2
    methods (Static)         
        %% > Wrap up Fig_2.
        function [] = WrapUp_Fig_2(Fig,inp,msh,blk)
            %  > Figure.
            figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
            Fig_2.Plot_1(inp,msh,blk,0);
        end
        
        %% > Auxiliary functions.
        % >> Plot 1.
        function [] = Plot_1(inp,msh,blk,len)
            %% > Cell borders.
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                patch(msh.c.XY_v{1,i}(:,1),msh.c.XY_v{1,i}(:,2),blk.f(i));
            end
            c = Fig_Tools.Colormap_cmocean();
            c.Label.String = '$\phi$';
            Fig_Tools.ChangeLook_1(inp,len);
        end
    end
end