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
            end
        end

        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1_1(msh,obj,x,fig,zoom)
            %  > Auxiliary variables.
            fid  = "V1_1_2D (1)";
            n    = x.a;
            L1   = Fig_V1_1_2D.Set_Legend_1;
          
            %  > Select variables.
            V{1} = msh.c.c.xy.v;
            V{2} = abs(obj.m{n}.At*(obj.x{1}.nv.x.c-obj.f.av.c));
            V{3} = msh.f.xy.v;            
            %  > Plot variables.
            Fig_Tools_2D.Var_3(1,V{1},V{2},1); 
            Fig_Tools_2D.Change_Colormap(fig,'thermal');
            
            Fig_Tools_2D.Var_1(1,V{3},"k","-",fig.LW);

            %  > Axis/legend,etc.
            Fig_Tools_2D.ChangeLook(fig,[V{3}{:}]');
            %  > Export(?)/Zoom(?).
            Fig_Tools_2D.Exp_Or_Zoom(fig,zoom,fid);
        end
        %  > 2.1.2. -------------------------------------------------------
        function [L] = Set_Legend_1()
            L{1} = join(["$|\bar{\tau}_{c}|$"]);
        end
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1_2(msh,obj,x,fig,zoom)
            %  > Auxiliary variables.
            fid  = "V1_1_2D (2)";
            n    = x.a;
            L1   = Fig_V1_1_2D.Set_Legend_1;
          
            %  > Select variables.
            V{1} = msh.c.c.xy.v;
            V{2} = abs(obj.x{1}.nv.x.c-obj.f.av.c);
            V{3} = msh.f.xy.v;            
            %  > Plot variables.
            Fig_Tools_2D.Var_3(1,V{1},V{2},1); 
            Fig_Tools_2D.Change_Colormap(fig,'thermal');
            
            Fig_Tools_2D.Var_1(1,V{3},"k","-",fig.LW);

            %  > Axis/legend,etc.
            Fig_Tools_2D.ChangeLook(fig,[V{3}{:}]');
            %  > Export(?)/Zoom(?).
            Fig_Tools_2D.Exp_Or_Zoom(fig,zoom,fid);
        end
    end
end