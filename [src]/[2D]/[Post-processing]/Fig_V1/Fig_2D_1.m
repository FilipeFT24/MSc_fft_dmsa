classdef Fig_2D_1
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,msh,obj)
            %  > Auxiliary variables.
            exp  = 0;
            run  = 0;
            zoom = 0;
            fig  = Fig_Tools.Set_fig(exp,run,zoom);

            if ~exp
                if inp.plot(1)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    %  > Diffusive term (x) - face "f".
                    f = 1;  %  >    f .
                    j = 2;  %  > (:,j).
                    k = 1;  %      {k}.
                    l = 1;  %      {l}.
                    Fig_2D_1.Plot_1(msh,f,obj.s{l}.sc{f,j}{k},obj.s{l}.sf{f,j}{k},fig);
                    subplot(1,2,2);
                    %  > Diffusive term (y) - face "f".
                    f = 2;  %  >    f .
                    j = 2;  %  > (:,j).
                    k = 2;  %      {k}.
                    l = 1;  %      {l}.
                    Fig_2D_1.Plot_1(msh,f,obj.s{l}.sc{f,j}{k},obj.s{l}.sf{f,j}{k},fig);
                end
            end
        end

        %% > 2. -----------------------------------------------------------
        function [] = Plot_1(msh,f,sc,sf,fig)
            %  > Auxiliary variables.
            fig.fid = "2D_1";
            LW      = 2.5;
            MS      = [3.5,7.0];
            
            %  > Axis/legend,etc.
            Fig_Tools.Map_2D_1(msh,fig);
            %  > Plot variables.
            Fig_Tools.Var_2D_1(msh.f.xy.v  (f),fig.C(1,:),"-",LW);
            Fig_Tools.Var_2D_2(msh.c.c.xy.c   ,"k"       ,"o",MS(1));
            for i = 1:numel(sc)
                Fig_Tools.Var_2D_2(msh.c.c.xy.c(sc{i},:),fig.C(i,:),"s",MS(2));
                Fig_Tools.Var_2D_3(msh.c.c.xy.v(sc{i})  ,fig.C(i,:)    ,fig.FA);
            end
            for i = 1:numel(sf)
                Fig_Tools.Var_2D_2(msh.f.xy.c  (sf{i},:),fig.C(i,:),"s",MS(2));
            end
            Fig_Tools.Export(fig);
        end
    end
end