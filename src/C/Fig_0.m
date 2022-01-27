classdef Fig_0
    methods (Static)
        %% > Wrap-up Fig_0.
        function [] = WrapUp_Fig_0(Plot_0,Exp_0,Fig,n)
            if Plot_0
                %  > Triangle (x,y) coordinates.
                [xy_T(1,:),xy_T(2,:),xy_T(3,:)] = ...
                    deal([-0.5,0.0],[0.5,0.0],[0.0,sqrt(3)./2]);                
                %  > Figure.
                figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
                subplot(1,2,1);
                Fig_0.Plot_1(n,xy_T);
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
        function [] = Plot_1(n,xy)
            %  > Quadrature abcissas/weights.
            Q = quadtriangle(n,'Domain',xy,...
                'Type','product','Symmetry','allowAsymmetric','Weights','allowNegative','Points','allowOutside');
            quadplot(Q,'PlotTitle','Off','PointLabels','Off');
        end
        % >> Plot 2.
        function [] = Plot_2(n)
            %  > Quadrature abcissas/weights.
            Q = quadsquare(n,'Type',...
                'productLegendre','Symmetry','allowAsymmetric','Weights','allowNegative','Points','allowOutside');
            quadplot(Q,'PlotTitle','Off','PointLabels','Off'); set(gcf,'color','w');
        end
    end
end