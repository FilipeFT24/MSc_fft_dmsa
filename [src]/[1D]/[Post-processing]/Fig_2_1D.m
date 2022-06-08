classdef Fig_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(flag,msh,obj)
            %  > Auxiliary variables.
            exp     = 0;
            run     = 0;
            zoom    = 0;
            fig     = Fig_Tools.Set_fig(exp,run,zoom);
            fig.dir = '[Post-processing]/[.pdf files]';
            
            %  > Plot(?).
            if flag
                if ~exp
                    %  > Select variables.
                    j    = 1;
                    A{1} = obj.e.a.t.f_abs;
                    A{2} = obj.e.a.n_abs.t.f(j,:);
                    B{1} = [obj.e.a.c.c_abs,obj.e.a.t.c_abs];
                    B{2} = [obj.e.a.n_abs.c(j),obj.e.a.n_abs.t.c(j)];
                    %  > Plot/export...
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    fig.fid = "1D_2 (1)";
                    Fig_2_1D.Plot_1(msh.f.Xv,msh.f.Xv,A,fig,[1,0,j]);
                    subplot(1,2,2);
                    fig.fid = "1D_2 (2)";
                    Fig_2_1D.Plot_1(msh.c.Xc,msh.f.Xv,B,fig,[2,0,j]);
                end
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1(X,Y,Z,fig,flag)
            %  > Auxiliary variables.
            L  = Fig_2_1D.Set_L(flag);
            M1 = repelem(":o",size(Z{1},2)); 
            M2 = repelem(" -",size(Z{1},2));
            %  > Plot variables.
            [L1,P1,Y1] = Fig_Tools.Var_1D_1(fig,M1,L(1,:),X,Z{1});
            [L2,P2,Y2] = Fig_Tools.Var_1D_3(fig,M2,L(2,:),X,Z{2});
            %  > Axis/legend,etc.
            Fig_Tools.Map_1D_1(fig,1,X,fig.L);
            Fig_Tools.Map_1D_2(fig,[L1,L2],[P1,P2],1,Y,[Y1;Y2],[-1,1.5],2)
            Fig_Tools.Export  (fig);
        end
        %  > 2.1.1. -------------------------------------------------------
        function [L] = Set_L(flag)
            S1(1) = Fig_Tools.Set_str_2(flag(3));
            switch flag(1)
                case 1
                    S2(1) = Fig_Tools.Set_str_1(1);
                    S2(2) = Fig_Tools.Set_str_1(2);
                    S2(3) = Fig_Tools.Set_str_1(3);
                    switch flag(2)
                        case false
                            %  > w/o predicted.
                            L{1,1} = join(["$ |\bar{\tau}_{f_{\left(a\right)}}^{",S2(1),"}|$"]);
                            L{1,2} = join(["$ |\bar{\tau}_{f_{\left(a\right)}}^{",S2(2),"}|$"]);
                            L{1,3} = join(["$ |\bar{\tau}_{f_{\left(a\right)}}^{",S2(3),"}|$"]);
                            L{2,1} = join(["$\|\bar{\tau}_{f_{\left(a\right)}}^{",S2(1),"}\|_{_{",S1(1),"}}$"]);
                            L{2,2} = join(["$\|\bar{\tau}_{f_{\left(a\right)}}^{",S2(2),"}\|_{_{",S1(1),"}}$"]);
                            L{2,3} = join(["$\|\bar{\tau}_{f_{\left(a\right)}}^{",S2(3),"}\|_{_{",S1(1),"}}$"]);
                        case true
                            %  > w/  predicted.
                        otherwise
                            return;
                    end
                case 2
                    switch flag(2)
                        case false
                            %  > w/o predicted.
                            L{1,1} = join(["$|e_{c_{\left(a\right)}}|$"]);
                            L{1,2} = join(["$|\bar{\tau}_{c_{\left(a\right)}}|$"]);
                            L{2,1} = join(["$\|e_{c_{\left(a\right)}}\|_{_{",S1(1),"}}$"]);
                            L{2,2} = join(["$\|\bar{\tau}_{c_{\left(a\right)}}\|_{_{",S1(1),"}}$"]);
                        case true
                            %  > w/  predicted.
                        otherwise
                            return;
                    end
                otherwise
                    return;
            end
        end
    end
end