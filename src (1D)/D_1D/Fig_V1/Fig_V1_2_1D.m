classdef Fig_V1_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,obj,msh)
            %  > Auxiliary variables.
            exp = 0;
            fig = Fig_Tools_1D.Set_fig(0,exp);
            m   = 1;
            n   = size(obj.e.t.f_abs,2);
            
            if ~exp
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                for i = 1:n
                    subplot(1,n,i);
                    str(i) = Fig_Tools_1D.Set_str_2(i);
                    Fig_V1_2_1D.Plot_1_1(obj.e.a{m}.t.f_abs(:,i),obj.e.t.f_abs{i},inp.ps.p(i),msh,fig,str(i),0);
                end
                %   figure; set(gcf,'Units','pixels','Position',fig.Position);
                %   Fig_V1_2_1D.Plot_1_2(obj.e.a{m}.t.c_abs,obj.e.t.c_abs,inp.ps.p(i),msh,fig,0);  
            else
                plot_1 = [0,1];
                if plot_1(1)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    str = Fig_Tools_1D.Set_str_2(1);
                    Fig_V1_2_1D.Plot_1_1(obj.e.a{m}.t.f_abs(:,1),obj.e.t.f_abs{1},inp.ps.p(1),msh,fig,str,0);
                end
                if plot_1(2)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    str = Fig_Tools_1D.Set_str_2(2);
                    Fig_V1_2_1D.Plot_1_1(obj.e.a{m}.t.f_abs(:,2),obj.e.t.f_abs{2},inp.ps.p(2),msh,fig,str,0);
                end
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1_1(et_a,et_x,p,msh,fig,str,zoom)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            m     = size(et_a,2);
            n     = size(et_x,2);
            
            %  > Plot variables.
            %  > #1.
            for j = 1:m+n
                if j == 1
                    M1(j) = "-o";
                else
                    M1(j) = ":v";
                end
            end
            L1{1} = join(["$|\bar{\tau}_{f^{\left(a\right)}}^{",str,"}|$"]);
            for j = 1:n
                L1{j+1} = join(["$|\bar{\tau}_{f^{\left(",num2str(p+j-1),"\right)}}^{",str,"}|$"]);
            end
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,msh.f.Xv,[et_a,et_x]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,L1,P1,msh.f.Xv,Y1,[-1.5,0],1); 
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,"tau_f");
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_1_2(et_a,et_x,p,msh,fig,zoom)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            m     = size(et_a,2);
            n     = size(et_x,2);
            
            %  > Plot variables.
            %  > #1.
            for j = 1:m+n
                if j == 1
                    M1(j) = "-o";
                else
                    M1(j) = ":v";
                end
            end
            L1{1} = "$|\bar{\tau}_{c^{\left(a\right)}}|$";
            for j = 1:n
                L1{j+1} = join(["$|\bar{\tau}_{c^{\left(",num2str(p+j-1),"\right)}}|$"]);
            end
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,msh.c.Xc,[et_a,et_x]);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,L1,P1,msh.c.Xc,Y1,[-1.5,0],1); 
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,"tau_c");
        end
    end
end