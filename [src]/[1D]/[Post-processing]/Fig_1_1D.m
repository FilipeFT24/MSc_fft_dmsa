classdef Fig_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(flag,msh)
            %  > Auxiliary variables.
            exp     = 0;
            run     = 0;
            zoom    = 0;
            fig     = Fig_Tools.Set_fig(exp,run,zoom);
            fig.fid = "1D_1";
            fig.dir = '[Post-processing]/[.pdf files]';
            
            %  > Plot(?).
            if flag
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                %  > Plot variables.
                hold on;
                plot(msh.f.Xv,repelem(0,msh.f.Nf),'-','Color','k','MarkerFaceColor','k','MarkerSize',15,'Linewidth',0.25);
                plot(msh.f.Xv,repelem(0,msh.f.Nf),'|','Color','k','MarkerFaceColor','k','MarkerSize',25,'Linewidth',0.25);
                %  > Axis/legend,etc.
                Fig_Tools.Map_1D_1(fig,1,msh.f.Xv,fig.L);
                axis off;
                Fig_Tools.Export  (fig);
            end
        end
    end
end