classdef Fig_V1_1_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(plot,msh,obj)
            %  > Auxiliary variables.
            exp = 0;
            fig = Fig_Tools_2D.Set_fig(0,exp);
            x.a = 1; %  > n.
            
            if ~exp
                if plot(1)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    Fig_V1_1_2D.Plot_1_1(msh,obj,x,fig,0);
                    subplot(1,2,2);
                    Fig_V1_1_2D.Plot_1_2(msh,obj,x,fig,0);
                end
            else
                if plot(1)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    Fig_V1_1_2D.Plot_1_1(msh,obj,x,fig,0);
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    Fig_V1_1_2D.Plot_1_2(msh,obj,x,fig,0);
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1_1(msh,obj,x,fig,zoom)
            %  > Auxiliary variables.
            fid             = "V1_1_2D (1)";
            n               = x.a;
            v               = obj.e.a{n}.t.c_abs;
            y.str           = 'thermal';
            y.title         = Fig_V1_1_2D.Set_Legend(1);
            [y.c(1),y.c(2)] = MinMaxElem(v);
            
            %  > Select variables.
            V{1} = msh.c.c.xy.v;
            V{2} = v;
            %  > Plot variables.
            Fig_Tools_2D.Var_3(V{1},V{2},1);
            Fig_Tools_2D.ChangeLook_2(fig,y);
            
            %  > Axis/legend,etc.
            Fig_Tools_2D.ChangeLook_1(fig,msh);
            %  > Export(?)/Zoom(?).
            Fig_Tools_2D.Exp_Or_Zoom (fig,zoom,fid);
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_1_2(msh,obj,x,fig,zoom)
            %  > Auxiliary variables.
            fid             = "V1_1_2D (2)";
            n               = x.a;
            v               = obj.e.a{n}.c.c_abs;
            y.str           = 'thermal';
            y.title         = Fig_V1_1_2D.Set_Legend(2);
            [y.c(1),y.c(2)] = MinMaxElem(v);
            
            %  > Select variables.
            V{1} = msh.c.c.xy.v;
            V{2} = v;
            %  > Plot variables.
            Fig_Tools_2D.Var_3(V{1},V{2},1);
            Fig_Tools_2D.ChangeLook_2(fig,y);
            
            %  > Axis/legend,etc.
            Fig_Tools_2D.ChangeLook_1(fig,msh);
            %  > Export(?)/Zoom(?).
            Fig_Tools_2D.Exp_Or_Zoom (fig,zoom,fid);
        end
        % >> 2.3 -------------------------------------------------------
        function [L] = Set_Legend(t)
            switch t
                case 1
                    L = join(["$|\bar{\tau}_{c^{\left(a\right)}}|$"]);
                case 2
                    L = join(["$|e_{c^{\left(a\right)}}|$"]);
                otherwise
                    return;
            end
        end
    end
end