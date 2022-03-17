classdef Fig_4_1D
    methods (Static)
        function [] = Plot(inp,msh,pde)
            %  > Auxiliary variables.
            Exp = 0;
            fig = Fig_Tools_1D.Set_fig(Exp);
            p   = A_2_1D.Compute_p(inp.ee.p,inp.ee.s);
            n   = 2;
            N   = [1,2];
            x   = [1,3];
            if length(x) ~= 2
                return;
            end
            
            if ~Exp
                plot = [0,1];
                if plot(1)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    for i = 1:n
                        subplot(1,n,i);
                        Fig_4_1D.Plot_1(msh,pde,i,p(x,i),x,fig,Exp,0);
                    end
                end
                if plot(2)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    for i = 1:n
                        subplot(1,n,i);
                        Fig_4_1D.Plot_2(msh,pde,i,p(x,i),x,fig,Exp,0);
                    end
                end
            else
                for i = 1:n
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    Fig_4_1D.Plot_1(msh,pde,i,p(x,i),x,fig,Exp,Zoom);
                end
            end
        end
        % >> 1. -----------------------------------------------------------
        function [] = Plot_1(msh,pde,i,p,x,fig,Exp,Zoom)
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C       = linspecer(9,'qualitative');
            str         = Fig_Tools_1D.Switch_Legend(i);
            k           = length(x);
            [l(1),l(2)] = MinMaxElem(p);
            
            for j = 1:k
                fig.M     (j) = "-v";
                fig.L     {j} = join(["$|\tau_{f}^{",str,"\left(",num2str(p(j)),"\right)}|$"]);
                Var_1   (:,j) = pde.et.av{j}(:,i);
                fig.M   (j+k) = "-.o";
                fig.L   {j+k} = join(["$|",str,"_{f\left(a-",num2str(p(j)),"\right)}^{\left(",num2str(l(1)),"\right)}-",str,"_{f\left(a-",num2str(p(j)),"\right)}^{\left(",num2str(l(2)),"\right)}|$"]);
                Var_1 (:,j+k) = pde.et.df.t{x(j),diff(x)}(:,i);
            end

            %  > Plot variables.
            [fig,P1,Y1] = Fig_Tools_1D.Var_1(fig,msh.f.Xv,Var_1);           
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,P1,Y1);
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(i)]),fig.Folder)
            end
            %  > Zoom(?).
            if Zoom
                zp = BaseZoom;
                zp.plot(Exp);
            end
        end
        % >> 2. -----------------------------------------------------------
        function [] = Plot_2(msh,pde,i,p,x,fig,Exp,Zoom)
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C       = linspecer(9,'qualitative');
            str         = Fig_Tools_1D.Switch_Legend(i);
            k           = length(x);
            [l(1),l(2)] = MinMaxElem(p);
            
            for j = 1:k+1
                if j ~= k+1
                    fig.M       (j) = "-v";
                    fig.L       {j} = join(["$|\tau_{f\left(a\right)}^{",str,"\left(",num2str(p(j)),"\right)}|$"]);
                    Var_1     (:,j) = pde.et.av{j}(:,i);
                    fig.M   (j+k+1) = "-.s";
                    fig.L   {j+k+1} = join(["$|\tau_{f}^{",str,"\left(",num2str(l(j)),"/",num2str(l(k+1-j)),"\right)}|$"]);
                    Var_1 (:,j+k+1) = pde.et.df.x{j,1}(:,i);
                else
                    fig.M       (j) = ":o";
                    fig.L       {j} = join(["$|\tau_{f\left(a\right)}^{",str,"\left(",num2str(l(1)),"\right)}-\tau_{f\left(a\right)}^{",str,"\left(",num2str(l(2)),"\right)}|$"]);
                    Var_1     (:,j) = pde.et.df.a{x(1),x(2)-1}(:,i);                    
                end
            end

            %  > Plot variables.
            [fig,P1,Y1] = Fig_Tools_1D.Var_1(fig,msh.f.Xv,Var_1);           
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,P1,Y1);
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(i)]),fig.Folder)
            end
            %  > Zoom(?).
            if Zoom
                zp = BaseZoom;
                zp.plot(Exp);
            end
        end
    end
end