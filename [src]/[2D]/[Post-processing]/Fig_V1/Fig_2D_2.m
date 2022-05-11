classdef Fig_2D_2
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,msh,obj)
            %  > Auxiliary variables.
            exp  = 0;
            run  = 0;
            zoom = 0;
            fig  = Fig_Tools.Set_fig(exp,run,zoom);
            n    = 1;
            
            if ~exp
                if inp.plot(2)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    Fig_2D_2.Plot_1(msh,obj.e.a{n}.t.c_abs,1,fig);
                    subplot(1,2,2);
                    Fig_2D_2.Plot_1(msh,obj.e.a{n}.c.c_abs,2,fig);
                end
            else
                if inp.plot(2)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    Fig_2D_2.Plot_1(msh,obj.e.a{n}.t.c_abs,1,fig);
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    Fig_2D_2.Plot_1(msh,obj.e.a{n}.c.c_abs,2,fig);
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1(msh,v,x,fig)
            %  > Auxiliary variables.
            fig.fid         = "2D_2";
            fig.LW          = 0.1;
            fig.FT_1        = 12.00;
            y.str           = 'thermal';
            y.title         = Fig_2D_2.Set_Legend(x);
            [y.c(1),y.c(2)] = MinMaxElem(v);
            
            %  > Select variables.
            V{1} = msh.c.c.xy.v;
            V{2} = v;
            %  > Plot variables.
            Fig_Tools.Var_2D_3(V{1},V{2},1);            
            %  > Axis/legend,etc.
            Fig_Tools.Map_2D_1(msh,fig);
            Fig_Tools.Map_2D_2(fig,y);
            Fig_Tools.Map_2D_3(fig);
        end
        % >> 2.2 -------------------------------------------------------
        function [L] = Set_Legend(t)
            switch t
                case 1, L = join(["$|\bar{\tau}_{c^{\left(a\right)}}|$"]);
                case 2, L = join(["$|e_{c^{\left(a\right)}}|$"]);
                otherwise
                    return;
            end
        end
    end
end