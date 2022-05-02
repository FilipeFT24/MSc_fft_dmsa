classdef Fig_V1_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(plot,inp,msh,obj)
            %  > Auxiliary variables.
            exp = 0;
            fig = Fig_Tools_1D.Set_fig(1,exp);
            x.a = 1; %  > (:,j): L_Norm.
            x.b = 1; %  >   {n}.
                       
            if ~exp
                if plot
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    n = 2;
                    for i   = 1:n
                        x.c = i;
                        subplot(1,n,i);
                        Fig_V1_2_1D.Plot_1_1(obj.e,obj.m,x,fig,0);
                    end
                end
            else
                if plot
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    x.c = 1;
                    Fig_V1_2_1D.Plot_1_1(obj.e,obj.m,x,fig,0);
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    x.c = 2;
                    Fig_V1_2_1D.Plot_1_1(obj.e,obj.m,x,fig,0);
                end
            end
        end

        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 2.1.1. -------------------------------------------------------
        function [] = Plot_1_1(obj_e,obj_m,x,fig,zoom)
            %  > Auxiliary variables.
            fid  = "V1_2_1D (1)";
            j    = x.a;
            n    = x.b;
            plot = [1,0];
            f    = ["a","da"];
            m    = length(f);
            L1   = Fig_V1_2_1D.Set_Legend_1(x.c,f(1),j);
            L2   = Fig_V1_2_1D.Set_Legend_1(x.c,f(2),j);
                        
            %  > Select variables.
            for i = 1:size(obj_m,1)
                nnz_all(i,1) = obj_m{i,1}.nnz.Ac(1);
                nnz_all(i,2) = obj_m{i,1}.nnz.Ac(2);
                nnz_all(i,3) = obj_m{i,1}.nnz.At;
            end
            switch x.c
                case 1
                    M1            = ["-o","-v","-d"];
                    M2            = [":o",":v",":d"];
                    nnz           = nnz_all;
                    [ys(1),ys(2)] = MinMaxElem(nnz);
                    for i = 1:size(obj_e,1)
                        for k = 1:size(nnz,2)
                            for l = 1:m
                                V{l}(i,k) = obj_e(i).(f(l)){n}.t.n_abs.f(j,k);
                            end
                        end
                    end
                case 2
                    M1            = ["-o","-v"];
                    M2            = [":o",":v"];
                    nnz           = [nnz_all(:,3),nnz_all(:,3)];
                    [ys(1),ys(2)] = MinMaxElem(nnz);
                    for i = 1:size(obj_e,1)
                        for k = 1:m
                            V{k}(i,1) = obj_e(i).(f(k)){n}.t.n_abs.c(j);
                            V{k}(i,2) = obj_e(i).(f(k)){n}.c.n_abs  (j);
                        end
                    end
                otherwise
                    return;
            end
            %  > Plot variables.
            [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,nnz,V{1});
            [L2,P2,Y2] = Fig_Tools_1D.Var_1(fig,M2,L2,nnz,V{2});
            
            %  > Axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,[L1,L2],[P1,P2],0,ys,[Y1;Y2],[-1,1],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,fid);
        end
        %  > 2.1.2. -------------------------------------------------------
        function [L] = Set_Legend_1(c,f,j)
            S(1) = Fig_Tools_1D.Set_str_1(1);
            S(2) = Fig_Tools_1D.Set_str_1(2);
            S(3) = Fig_Tools_1D.Set_str_3(j);
            switch c
                case 1
                    L{1} = join(["$\|\bar{\tau}_{f^{\left(",f,"\right)}}^{",S(1),"}\|_{_{",S(3),"}}$"]);
                    L{2} = join(["$\|\bar{\tau}_{f^{\left(",f,"\right)}}^{",S(2),"}\|_{_{",S(3),"}}$"]);
                    L{3} = join(["$\|\bar{\tau}_{f^{\left(",f,"\right)}}^{\phantom{\nabla\phi}}\|_{_{",S(3),"}}$"]);
                case 2
                    L{1} = join(["$\|\bar{\tau}_{c^{\left(",f,"\right)}}^{\phantom{\nabla\phi}}\|_{_{",S(3),"}}$"]);
                    L{2} = join(["$\|e_{c^{\left(",f,"\right)}}^{\phantom{\nabla\phi}}\|_{_{",S(3),"}}$"]);
                otherwise
                    return;
            end
        end
    end
end