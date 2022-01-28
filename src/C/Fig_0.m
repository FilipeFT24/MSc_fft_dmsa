classdef Fig_0
    methods (Static)
        %% > Wrap-up Fig_0.
        function [] = WrapUp_Fig_0(Plot_0,Exp_0,Fig,n)
            if Plot_0               
                %  > Figure.
                figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
                subplot(1,2,1);
                Fig_0.Plot_1(n);
                subplot(1,2,2);
                Fig_0.Plot_2(n);
                %  > Export as .pdf.
                if Exp_0
                    Fig_Tools.Export_PDF('Fig_0','../[Figures]/Fig_0');
                end
            end
        end
        
        %% > Auxiliary functions.
        % >> Plot 1.
        function [] = Plot_1(n)
            %  > Quadrature abcissas/weights.
            Q = quadtriangle(n,'Domain',[0,0;1,0;0,1]);
            quadplot(Q,'PlotTitle','Off','PointLabels','Off');
        end
        % >> Plot 2.
        function [] = Plot_2(n)
            %  > Quadrature abcissas/weights.
            Q = quadsquare(n);
            quadplot(Q,'PlotTitle','Off','PointLabels','Off'); set(gcf,'color','w');
        end
    end
end