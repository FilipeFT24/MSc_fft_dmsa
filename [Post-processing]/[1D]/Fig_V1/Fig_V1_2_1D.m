classdef Fig_V1_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,msh,obj,plot)
            %  > Auxiliary variables.
            exp = 0;
            fig = Fig_Tools_1D.Set_fig(0,exp);
            x.a = 1; %  > (:,j).
            x.b = 1; %  >   {n}.
            
            if ~exp
                if plot
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    x.a = 1;
                    Fig_V1_2_1D.Plot_1X(msh.f.Xv,obj,x,fig,0);
                    subplot(1,2,2);
                    x.a = 2;
                    Fig_V1_2_1D.Plot_1X(msh.f.Xv,obj,x,fig,0);
                end
            else
                if plot
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    x.a = 1;
                    Fig_V1_2_1D.Plot_1X(msh.f.Xv,obj,x,fig,0);
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    x.a = 2;
                    Fig_V1_2_1D.Plot_1X(msh.f.Xv,obj,x,fig,0);
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 2.1.1. -------------------------------------------------------
        function [] = Plot_1X(Xv,obj,x,fig,zoom)
            %  > Auxiliary variables.
            fid = "V1_2_1D (1)";
            j   = x.a;
            n   = x.b;
            
            %  > Select variables.
            L1      = Fig_V1_2_1D.Set_Legend_1(j);
            M1      = ["--o","--s","--v","--d","-.o","-.s","-.v"];
            f1      = ["a","da","p","d"];
            l1      = length(f1);
            for i = 1:l1
                V1(:,i) = obj.e.(f1(i)){n}.t.f_abs(:,j);         % \tau_f.
            end
            V1(:,5) = obj.x{n}  .xf.a(:,j)-obj.x{n}  .xf.x(:,j); % \X(n).
            V1(:,6) = obj.x{n+1}.xf.a(:,j)-obj.x{n+1}.xf.x(:,j); % \X(n+1).
            V1(:,7) = V1(:,5)-V1(:,6);                           % \X(dn). 
            V1      = abs(V1);
            %  > Plot variables.
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,Xv,V1);
            
            %  > Axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,L1,P1,1,Xv,Y1,[-1,1],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,fid);
        end
        %  > 2.1.2. -------------------------------------------------------
        function [L] = Set_Legend_1(j)
            S    = Fig_Tools_1D.Set_str_2(j);
            L{1} = join(["$|\bar{\tau}_{f^{\left(a\right)\phantom{d-p}}}^{",S,"}|$"]);
            L{2} = join(["$|\bar{\tau}_{f^{\left(da\right)\phantom{-p}}}^{",S,"}|$"]);
            L{3} = join(["$|\bar{\tau}_{f^{\left(p\right)\phantom{d-a}}}^{",S,"}|$"]);
            L{4} = join(["$|\bar{\tau}_{f^{\left(da-p\right)}}^{",S,"}|$"]);
            L{5} = join(["$|",S,"_{f_{\left(a-m\right)}}^{\left(m\right)\phantom{-n}}|$"]);
            L{6} = join(["$|",S,"_{f_{\left(a-m\right)}}^{\left(n\right)\phantom{-m}}|$"]);
            L{7} = join(["$|",S,"_{f_{\left(a-m\right)}}^{\left(m-n\right)}|$"]);           
        end
    end
end