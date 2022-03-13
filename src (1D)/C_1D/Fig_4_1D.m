classdef Fig_4_1D
    methods (Static)
        function [] = Plot(msh,pde,p)
            %  > Auxiliary variables.
            Exp = 0;
            fig = Fig_Tools_1D.Set_fig(Exp);
            n   = 2;
            N   = [1,2];
            x   = [1,2];
            if length(x) ~= 2
                return;
            end
            
            if ~Exp
                figure(N(1)); set(gcf,'Units','pixels','Position',fig.Position);
                for i = 1:n
                    subplot(1,n,i);
                    Fig_4_1D.Plot_1(msh,pde,i,p(x,i),x,fig,Exp);
                end
            else
                for i = 1:n
                    figure(N(i)); set(gcf,'Units','pixels','Position',fig.Position);
                    Fig_4_1D.Plot_1(msh,pde,i,p(x,i),x,fig,Exp);
                end
            end
        end
        % >> 1. -----------------------------------------------------------
        function [] = Plot_1(msh,pde,i,p,x,fig,Exp)
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C       = linspecer(9,'qualitative');
            str         = Fig_Tools_1D.Switch_Legend(i);
            k           = length(x);
            [l(1),l(2)] = MinMaxElem(p);
            
            for j = 1:k+1
                if j ~= k+1
                    %  > [1,2].
                    fig.M       (j) = "-v";
                    fig.L1      {j} = join(["$|\tau_{f}^{",str,"\left(",num2str(p(j)),"\right)}|$"]);
                    Var_1     (:,j) = pde.et.av{j}(:,i);
                    %  > [4,5].
                    fig.M   (j+k+1) = "-.s";
                    fig.L1  {j+k+1} = join(["$|",str,"_{f\left(a-",num2str(p(j)),"\right)}^{\left(",num2str(l(1)),"\right)}-",str,"_{f\left(a-",num2str(p(j)),"\right)}^{\left(",num2str(l(2)),"\right)}|$"]);
                    Var_1 (:,j+k+1) = pde.et.df.t{x(j),diff(x)}(:,i);
                else
                    %  > [3].
                    fig.M       (j) = ":o";
                    fig.L1      {j} = join(["$|\tau_{f}^{",str,"\left(",num2str(l(1)),"\right)}-\tau_{f}^{",str,"\left(",num2str(l(2)),"\right)}|$"]);
                    Var_1     (:,j) = pde.et.df.a{x(1),x(2)}(:,i);                    
                end
            end
            fig.L2{1} = "$x$";
            fig.L2{2} = "$\textrm{Error magnitude}$";

            %  > Plot variables.
            [fig,P1,Y1] = Fig_Tools_1D.Var_1(fig,msh.f.Xv,Var_1);           
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,P1,Y1);
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(i)]),fig.Folder)
            end
        end
    end
end