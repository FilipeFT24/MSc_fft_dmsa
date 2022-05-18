classdef Fig_2D_1
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [] = Plot(inp,msh,obj)
            %  > Auxiliary variables.
            exp  = 0;
            run  = 0;
            zoom = 0;
            fig  = Fig_Tools.Set_fig(exp,run,zoom);
            x.c  = 1; %  > n.
            
            if ~exp
                if inp.plot(1)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    x.a = 3;                             %  >    f .
                    x.b = 1;                             %  > (:,j).
                    Fig_2D_1.Plot_1(msh,obj,x,fig);
                    subplot(1,2,2);
                    x.a = 58;                             %  >    f .
                    x.b = 2;                             %  > (:,j).
                    Fig_2D_1.Plot_1(msh,obj,x,fig);
                end
            else
                if inp.plot(1)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    x.a = 1;                             %  >    f .
                    x.b = 1;                             %  > (:,j).
                    Fig_2D_1.Plot_1(msh,obj,x,fig);
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    x.a = 2;                             %  >    f .
                    x.b = 2;                             %  > (:,j).
                    Fig_2D_1.Plot_1(msh,obj,x,fig);
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        function [f] = Set_f(msh,ft)
            switch ft
                case "bnd"
                    v = find(~msh.f.logical);
                case "blk"
                    v = find( msh.f.logical);
                otherwise
                    return;
            end
            f = v(randperm(numel(v),1));
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1(msh,obj,x,fig)
            %  > Auxiliary variables.
            fid = "2D_1";
            f   = x.a;
            j   = x.b;
            n   = x.c;
            FA  = fig.FA;
            LW  = [2.5,0.1];
            MS  = [3.0,7.0];
            
            %  > Axis/legend,etc.
            Fig_Tools.Map_2D_1(msh,fig);
            %  > Select variables.
            V{1} = msh.f.xy.v;
            V{2} = msh.c.c.xy.c;
            for i = 1:numel(obj.s{n}.sc{f,j}{n})
                V{3}{i} = msh.c.c.xy.c(obj.s{n}.sc{f,j}{n}{i},:);
                V{4}{i} = msh.c.c.xy.v(obj.s{n}.sc{f,j}{n}{i});
            end
            for i = 1:numel(obj.s{n}.sf{f,j}{n})
                V{5}{i} = msh.f.xy.c  (obj.s{n}.sf{f,j}{n}{i},:);
            end
            %  > Plot variables.
            Fig_Tools.Var_2D_1(V{1}   ,"k"       ,"-",LW(2));
            Fig_Tools.Var_2D_2(V{2}   ,"k"       ,"o",MS(1));
            Fig_Tools.Var_2D_1(V{1}(f),fig.C(1,:),"-",LW(1));
            for i = 1:numel(obj.s{n}.sc{f,j}{n})
                Fig_Tools.Var_2D_2(V{3}{i},fig.C(i,:),"s",MS(2));
                Fig_Tools.Var_2D_3(V{4}{i},fig.C(i,:)    ,FA);
            end
            for i = 1:numel(obj.s{n}.sf{f,j}{n})
                Fig_Tools.Var_2D_2(V{5}{i},fig.C(i,:),"s",MS(2));
            end
            Fig_Tools.Map_2D_3(fig);
        end
    end
end