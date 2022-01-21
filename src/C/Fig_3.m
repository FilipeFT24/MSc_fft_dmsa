classdef Fig_3
    methods (Static)
        %% > Wrap-up Fig_3.
        function [] = WrapUp_Fig_3(Plot_3,Exp_3,Fig,inp,msh,pde,len)
            if Plot_3
                %  > Figure.
                figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
                Fig_3.Plot_1(inp,msh,pde,len);
                %  > Export as .pdf.
                if Exp_3
                    Fig_Tools.Export_PDF('Fig_3','../[Figures]/Fig_3');
                end
            end
        end
        
        %% > Auxiliary functions.
        % >> Plot 1.
        function [] = Plot_1(inp,msh,pde,len)           
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                patch(msh.c.xy_v{i}(:,1),msh.c.xy_v{i}(:,2),pde.c.F_Vol(i));
            end
            c = Fig_Tools.Colormap_cmocean('thermal');
            c.Label.String = '$S^\phi$';
            Fig_Tools.ChangeLook_1(inp,len);
            AdvancedColormap('thermal');           
        end
    end
end