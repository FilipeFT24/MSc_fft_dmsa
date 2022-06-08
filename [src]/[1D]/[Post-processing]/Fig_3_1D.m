classdef Fig_3_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(flag,msh,obj)
            %  > Auxiliary variables.
            exp     = 0;
            run     = 1;
            zoom    = 0;
            fig     = Fig_Tools.Set_fig(exp,run,zoom);
            fig.dir = '[Post-processing]/[.pdf files]';
            
            %  > Plot(?).
            if flag
                %  > Select variables.
                for i = 1:size(obj,1)
                    for j = 1:size(obj,2)
                        V1{1}(i,j) = obj(i,j).e.a.n_abs.t.f(1,3);
                        V1{2}(i,j) = obj(i,j).e.a.n_abs.c  (1);
                        V3{1}(i,j) = obj(i,j).e.a.n_abs.t.f(3,3);
                        V3{2}(i,j) = obj(i,j).e.a.n_abs.c  (3);
                    end
                    NNZ(i,1) = obj(i,1).m{1}.nnz.At;
                    h  (i,1) = msh(i,1).d.h;
                end
                if ~exp
                    %  > Plot...
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    Fig_3_1D.Plot_1(h,NNZ,V1,fig,[0,1]);
                    subplot(1,2,2);
                    Fig_3_1D.Plot_1(h,NNZ,V3,fig,[0,3]);
                else
                    %  > Plot/export...
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    fig.fid = "1D_2 (1)";
                    Fig_3_1D.Plot_1(h,NNZ,V1,fig,[0,1]);
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    fig.fid = "1D_2 (2)";
                    Fig_3_1D.Plot_1(h,NNZ,V3,fig,[0,3]);
                end
            end
        end
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1(h,X,Y,fig,flag)
            %  > Auxiliary variables.
            L = Fig_3_1D.Set_L(flag);
            M = ["-o","-o","-o";":v",":v",":v"];
            for i = 1:numel(Y)
                [LP{i},PP{i},YP{i}] = Fig_Tools.Var_1D_1(fig,M(i,:),L(i,:),X,Y{i});
            end
            %  > Axis/legend,etc.
            Fig_Tools.Map_1D_2(fig,[LP{:}],[PP{:}],0,X,[YP{:}],[-0.5,1.5],2)
            Fig_Tools.Export  (fig);
            %  > Slope.
            S1 = src_Tools.Slope(h,Y{1});
            S2 = src_Tools.Slope(h,Y{2});
        end
        %  > 2.1.1. -------------------------------------------------------
        function [L] = Set_L(flag)
            S1(1) = Fig_Tools.Set_str_2(flag(2));
            switch flag(1)
                case 0
                    %  > w/o predicted.
                    f = ["a","b","c"];
                    for i = 1:numel(f)
                        L{1,i} = join(["$\|\bar{\tau}_{f_{\left(",f(i),"\right)}}^{\phi}\|_{_{",S1(1),"}}$"]);
                        L{2,i} = join(["$\|e_{c_{\left(",f(i),"\right)}}\|_{_{",S1(1),"}}$"]);
                    end
                case 1
                    %  > w/  predicted.
                otherwise
                    return; 
            end
        end
    end
end