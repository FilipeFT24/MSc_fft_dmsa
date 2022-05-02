classdef Fig_V1_0_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [] = Plot(plot,msh,obj,f_type)
            %  > Auxiliary variables.
            exp = 0;
            fig = Fig_Tools_2D.Set_fig(0,exp);
            x.c = 1; %  > n.
            
            if ~exp
                if plot(1)
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    x.a = 1;%Fig_V1_0_2D.Set_f(msh,f_type); %  >    f .
                    x.b = 1;                             %  > (:,j).
                    Fig_V1_0_2D.Plot_1(msh,obj,x,fig,0);
                    subplot(1,2,2);
                    x.a = Fig_V1_0_2D.Set_f(msh,f_type); %  >    f .
                    x.b = 2;                             %  > (:,j).
                    Fig_V1_0_2D.Plot_1(msh,obj,x,fig,0);
                end
            else
            end
        end
        % >> 1.2. ---------------------------------------------------------
        function [f] = Set_f(msh,str)
            switch str
                case "bnd"
                    bnd = find(~msh.f.logical);
                    f   = bnd ( randperm(numel(bnd),1));
                case "blk"
                    blk = find( msh.f.logical);
                    f   = blk ( randperm(numel(blk),1));
                otherwise
                    return;
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1(msh,obj,x,fig,zoom)
            %  > Auxiliary variables.
            fid  = "V1_1_2D (1/2)";
            f    = x.a;
            j    = x.b;
            n    = x.c;
            
            %  > Select variables.
            V{1} = msh.f.xy.v;
            V{2} = msh.c.c.xy.c;
            V{3} = obj.s{n}.xt{f,j};
            %  > Plot variables.
            Fig_Tools_2D.Var_1(1,V{1}   ,"k","-",fig.LW);
            Fig_Tools_2D.Var_2(1,V{2}   ,"k","o",fig.MS./2.0);
            Fig_Tools_2D.Var_1(1,V{1}(f),fig.C(1,:),"-",fig.LW.*25);
            Fig_Tools_2D.Var_2(1,V{3}   ,fig.C(1,:),"s",fig.MS);
            
            %  > Axis/legend,etc.
            Fig_Tools_2D.ChangeLook (fig,[msh.f.xy.v{:}]');
            %  > Export(?)/Zoom(?).
            Fig_Tools_2D.Exp_Or_Zoom(fig,zoom,fid);
        end
    end
end