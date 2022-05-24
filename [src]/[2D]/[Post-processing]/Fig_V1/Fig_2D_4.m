classdef Fig_2D_4
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [] = Plot(inp,msh,obj)
            %  > Auxiliary variables.
            exp  = 0;
            run  = 0;
            zoom = 0;
            fig  = Fig_Tools.Set_fig(exp,run,zoom);
            
            if ~exp
                if inp.plot(3)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    f = 10; %  >    f .
                    j = 2;  %  > (:,j).
                    k = 1;  %      {k}.
                    l = 1;  %      {l}.
                    Fig_2D_4.Plot_1(msh,f,obj.e.f{f,j}{k},obj.s{l}.sc{f,j}{k},obj.s{l}.sf{f,j}{k},obj.s{l}.logical{f,j}{k},fig);
                    subplot(1,2,2);
                    f = 11; %  >    f .
                    j = 2;  %  > (:,j).
                    k = 2;  %      {k}.
                    l = 1;  %      {l}.
                    Fig_2D_4.Plot_1(msh,f,obj.e.f{f,j}{k},obj.s{l}.sc{f,j}{k},obj.s{l}.sc{f,j}{k},obj.s{l}.logical{f,j}{k},fig);
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
        function [] = Plot_1(msh,f,obj_e,sc,sf,logical,fig)
            %  > Auxiliary variables.
            fig.fid         = "2D_3";
            fig.LW          = 0.1;
            fig.FT_1        = 12.00;
            y.str           = 'thermal';
            y.title         = 'X';
            [y.c(1),y.c(2)] = MinMaxElem(obj_e);
            
            %  > Plot variables.
            Fig_Tools.Var_2D_3(msh.c.c.xy.v(cat(1,sc{:})),obj_e(logical),1);            
            %  > Axis/legend,etc.
            Fig_Tools.Map_2D_1(msh,fig);
            Fig_Tools.Map_2D_2(fig,y);
            Fig_Tools.Map_2D_3(fig);
            Fig_Tools.Var_2D_1(msh.f.xy.v(f),fig.C(1,:),"-",2.5);
        end
    end
end