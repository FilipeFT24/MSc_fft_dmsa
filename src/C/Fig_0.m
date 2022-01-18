classdef Fig_0
    methods (Static)
        %% > Wrap up Fig_0.
        function [] = WrapUp_Fig_0(Fig,n)
            %  > Triangle (x,y) coordinates.
            [xy_T(1,:),xy_T(2,:),xy_T(3,:)] = ...
                deal([-0.5,0.0],[0.5,0.0],[0.0,sqrt(3)./2]);
            
            %  > Figure.
            figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
            Fig_0.Plot_1(n,xy_T);
            Fig_0.Plot_2(n);
        end
        
        %% > Auxiliary functions.
        % >> Plot 1.
        function [] = Plot_1(n,xy)
            subplot(1,2,1);
            %  > Quadrature abcissas/weights.
            Q = quadtriangle(n,'Domain',xy,...
                'Type','product','Symmetry','allowAsymmetric','Weights','allowNegative','Points','allowOutside');
            %  > Plot...
            quadplot(Q,'PlotTitle','Off','PointLabels','Off')
        end
        % >> Plot 2.
        function [] = Plot_2(n)
            subplot(1,2,2);
            %  > Quadrature abcissas/weights.
            Q = quadsquare(n,...
                'Type','productLegendre','Symmetry','allowAsymmetric','Weights','allowNegative','Points','allowOutside');
            %  > Plot...
            quadplot(Q,'PlotTitle','Off','PointLabels','Off')
            set(gcf,'color','w');
        end
    end
end