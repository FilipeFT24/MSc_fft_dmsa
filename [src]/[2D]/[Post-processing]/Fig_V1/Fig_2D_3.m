classdef Fig_2D_3
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(msh,obj)
            %  > Auxiliary variables.
            exp  = 0;
            run  = 1;
            zoom = 0;
            fig  = Fig_Tools.Set_fig(exp,run,zoom);
            x.a  = 1; %  > (:,j).
            x.b  = 1; %  >   {n}.
            
            if ~exp
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                subplot(1,2,1);
                Fig_2D_3.Plot_1(msh,obj,x,fig);
                subplot(1,2,2);
                Fig_2D_3.Plot_2(msh,obj,x,fig);
            else
                figure; set(gcf,'Units','pixels','Position',fig.Position);
                Fig_2D_3.Plot_2(msh,obj,x,fig);
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 2.1.1. -------------------------------------------------------
        function [] = Plot_1(msh,obj,x,fig)
            %  > Auxiliary variables.
            fig.fid = "2D_3 (1)";
            j       = x.a;
            n       = x.b;
            
            %  > Select variables.
            for i = 1:size(obj,2)
                h   (i,1) = msh(i).d.h;
                for k = 1:numel(obj(i).m{n}.nnz.Ac)+1
                    if k ~= numel(obj(i).m{n}.nnz.Ac)+1
                        NNZ(i,k) = obj(i).m{n}.nnz.Ac(k);
                    else
                        NNZ(i,k) = obj(i).m{n}.nnz.At;
                    end
                end
                for k = ["x","y"]
                    V{1}.(k)(i,:) = obj(i).e.a{n}.n_abs.t.f.(k)(j,:);
                end
            end
            %  > Plot variables.
            L1         = Fig_2D_3.Set_Legend_1(j,1);
            M1         = ["--o",":o","-.o"];
            [L1,P1,Y1] = Fig_Tools.Var_1D_1(fig,M1,L1,NNZ,V{1}.x);
            L2         = Fig_2D_3.Set_Legend_1(j,2);
            M2         = ["--v",":v","-.v"];
            [L2,P2,Y2] = Fig_Tools.Var_1D_1(fig,M2,L2,NNZ,V{1}.y);
            %  > Axis/legend,etc.
            Fig_Tools.Map_1D_2(fig,[L1,L2],[P1,P2],0,NNZ,[Y1;Y2],[-1,0],2);
            
            
            Tools_1.Slope(h,V{1}.x)
            Tools_1.Slope(h,V{1}.y)
            
            
        end
        %  > 2.1.2. -------------------------------------------------------
        function [L] = Set_Legend_1(j,n)
            S(1) = Fig_Tools.Set_str_1(1);
            S(2) = Fig_Tools.Set_str_1(2);
            S(3) = Fig_Tools.Set_str_1(3);
            S(4) = Fig_Tools.Set_str_3(j);
            switch n
                case 1, d = "x";
                case 2, d = "y";
                otherwise
                    return;
            end
            L{1} = join(["$\|\bar{\tau}_{f^{\left(a\right)}_{",d,"}}^{",S(1),"}\|_{_{",S(4),"}}$"]);
            L{2} = join(["$\|\bar{\tau}_{f^{\left(a\right)}_{",d,"}}^{",S(2),"}\|_{_{",S(4),"}}$"]);
            L{3} = join(["$\|\bar{\tau}_{f^{\left(a\right)}_{",d,"}}^{",S(3),"}\|_{_{",S(4),"}}$"]);
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_2(msh,obj,x,fig)
            %  > Auxiliary variables.
            fig.fid = "2D_3 (2)";
            j       = x.a;
            n       = x.b;
            
            %  > Select variables.
            for i = 1:size(obj,2)
                h     (i,1) = msh(i).d.h;
                NNZ   (i,1) = obj(i).m{n}.nnz.At;
                V  {1}(i,1) = obj(i).e.a{n}.n_abs.c  (j);
                V  {2}(i,1) = obj(i).e.a{n}.n_abs.t.c(j);
            end
            %  > Plot variables.
            L1         = Fig_2D_3.Set_Legend_2(j);
            M1         = ["--o",":v"];
            [L1,P1,Y1] = Fig_Tools.Var_1D_1(fig,M1,L1,NNZ,[V{1},V{2}]);
            %  > Axis/legend,etc.
            Fig_Tools.Map_1D_2(fig,L1,P1,0,NNZ,Y1,[-1,0],2);
            
            
            Tools_1.Slope(h,V{1})
            Tools_1.Slope(h,V{2})
            
        end
        % >> 2.3 ----------------------------------------------------------
        function [L] = Set_Legend_2(j)
            S(1) = Fig_Tools.Set_str_3(j);
            L{1} = join(["$\|e_{c^{\left(a\right)}}\|_{_{",S(1),"}}$"]);
            L{2} = join(["$\|\bar{\tau}_{c^{\left(a\right)}}\|_{_{",S(1),"}}$"]);
        end
    end
end