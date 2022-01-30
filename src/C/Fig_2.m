classdef Fig_2
    methods (Static)
        %% > Wrap up-Fig_2.
        function [] = WrapUp_Fig_2(Plot_2,Exp_2,Fig,inp,msh,pde,len)
            if Plot_2
                %  > Figure.
                figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
                Fig_2.Plot(inp,msh,pde,len,Exp_2);
                %  > Export as .pdf.
                if Exp_2
                    Fig_Tools.Export_PDF('Fig_2','../[Figures]/Fig_2');
                end
            end
        end
        
        %% > Auxiliary functions.
        function [] = Plot(inp,msh,pde,len,Flag)
            %% > Cell borders.
            %  > NOTE: Add "'Linestyle','None'" to remove cell border.
            hold on;
            for i = 1:msh.c.NC
                patch(msh.c.xy_v{1,i}(:,1),msh.c.xy_v{1,i}(:,2),pde.E.EA(i));
            end
            c = Fig_Tools.Colormap_cmocean('thermal');
            c.Label.String = '$\textrm{Absolute error}\left(\epsilon_{abs}\right)$';
            Fig_Tools.ChangeLook_1(inp,len,Flag);
            AdvancedColormap('thermal');
        end
    end
end