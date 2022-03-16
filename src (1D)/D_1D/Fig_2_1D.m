classdef Fig_2_1D
    methods (Static)
        function [] = Plot(msh,pde,flag)
            if flag
                %  > Auxiliary variables.
                Exp = 0;
                fig = Fig_Tools_1D.Set_fig(Exp);
                n   = size(pde.e.t.a,2);
                N   = 1:n;
                
                if ~Exp
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    for i = 1:n
                        subplot(1,n,i);
                        Fig_2_1D.Plot_1(msh,pde,i,fig,Exp);
                    end
                else
                    for i = 1:n
                        figure; set(gcf,'Units','pixels','Position',fig.Position);
                        Fig_2_1D.Plot_1(msh,pde,i,fig,Exp);
                    end
                end
            end
        end
        % >> 1. -----------------------------------------------------------
        function [] = Plot_1(msh,pde,i,fig,Exp)
            %  > Auxiliary variables (colors/labels/etc.).
            fig.C = linspecer(9,'qualitative');
            str   = Fig_Tools_1D.Switch_Legend(i);
            n     = size(pde.av.f,2);
            for j = 1:size(pde.e.t.a{i},2)+1
                if j == 1
                    fig.M (j) = "-v";
                    fig.L1{j} = join(["$|\bar{\tau}_{f}^{",str,"\phantom{\left(j\right)}}|$"]);
                else
                    fig.M (j) = ":o";
                    fig.L1{j} = join(["$|\bar{\tau}_{f}^{",str,"\left(",num2str(msh.s.stl.p(i*n-1)-i+j),"\right)}|$"]);
                end
            end
            fig.L2{1} = "$x$";
            fig.L2{2} = "$\textrm{Error magnitude}$";
            
            %  > Plot variables.
            for j = 1:size(pde.e.t.a{i},2)+1
                if j == 1
                    Var_1(:,j) = pde.e.t.f_abs(:,i);
                else
                    Var_1(:,j) = pde.e.t.a{i}(:,j-1);
                end
            end
            [fig,P,Y] = Fig_Tools_1D.Var_1(fig,msh.f.Xv,Var_1);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,P,Y);
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF(join(["Fig_1_",num2str(i)]),fig.Folder)
            end
        end
    end
end