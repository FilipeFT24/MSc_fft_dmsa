classdef Fig_V1_0_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(msh)
            %  > Auxiliary variables.
            exp  = 0;
            fig  = Fig_Tools_1D.Set_fig(0,exp);
            
            %  > #1.
            figure; set(gcf,'Units','pixels','Position',fig.Position);
            Fig_V1_0_1D.Plot_1(msh,fig);
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1(msh,fig)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            if ~fig.exp
                a = 3.5;
                b = 1;
            else
                a = 5.0;
                b = 1;
            end
            
            %  > Plot variables.
            hold on;
            plot(msh.f.Xv,repelem(-b,msh.f.Nf),'-','Color','k','MarkerFaceColor','k','MarkerSize',   fig.MS,'Linewidth',fig.LW./25.0);
            plot(msh.f.Xv,repelem(-b,msh.f.Nf),'|','Color','k','MarkerFaceColor','k','MarkerSize',a.*fig.MS,'Linewidth',fig.LW./25.0);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot_1(fig,b,msh.f.Xv)
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,0,"Fig_V1_0_1");
        end
    end
end