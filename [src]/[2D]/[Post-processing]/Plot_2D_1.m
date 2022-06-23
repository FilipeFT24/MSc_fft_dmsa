classdef Plot_2D_1
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,msh,obj)
            %  > Auxiliary variables.
            x.e = 0; %  > Export.
            x.l = 3; %  > Label.
            x.z = 0; %  > Zoom.
            fig = Fig_Tools.Set_fig(0,x);
            
            %  > Select variables.
            if inp.plot{1}(1)
                %  > Diffusive term (x) - face "f".
                f(1) = 22; %  >    f .
                j(1) = 2; %  > (:,j).
                k(1) = 1; %      {k}.
                l(1) = 1; %      {l}.
                %  > Diffusive term (y) - face "f".
                f(2) = 23; %  >    f .
                j(2) = 2; %  > (:,j).
                k(2) = 2; %      {k}.
                l(2) = 1; %      {l}.
            end
            if inp.plot{1}(2)
                V = [obj.e.a.t.c_abs,obj.e.a.c.c_abs];
                y = Plot_2D_1.Set_y(V,1);
            end
            %  > Plot variables.
            if inp.plot{1}(1)
                figure; set(gcf,'Units','pixels','Position',fig.pos);
                subplot(1,2,1);
                Plot_2D_1.Plot_1(fig,inp,msh,...
                    f(1),obj.s{l(1)}.u.p{j(1)}{k(1)}(f(1),:),obj.s{l(1)}.sc{f(1),j(1)}{k(1)},obj.s{l(1)}.sf{f(1),j(1)}{k(1)},obj.f.qd.xu);
                subplot(1,2,2);
                Plot_2D_1.Plot_1(fig,inp,msh,...
                    f(2),obj.s{l(2)}.u.p{j(2)}{k(2)}(f(2),:),obj.s{l(2)}.sc{f(2),j(2)}{k(2)},obj.s{l(2)}.sf{f(2),j(2)}{k(2)},obj.f.qd.xu);
                if x.e
                    Fig_Tools.SubPlot_pdf(fig.dir,'Plot_2D_1(1).pdf');
                end
            end
            if inp.plot{1}(2)
                figure; set(gcf,'Units','pixels','Position',fig.pos);
                subplot(1,2,1);
                Plot_2D_1.Plot_2(fig,msh,V(:,1),y{1});
                subplot(1,2,2);
                Plot_2D_1.Plot_2(fig,msh,V(:,2),y{2});
                if x.e
                    Fig_Tools.SubPlot_pdf(fig.dir,'Plot_2D_1(2).pdf');
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1(fig,inp,msh,f,p,sc,sf,xu)
            %  > Auxiliary variables.
            lw     = 2.5;
            ms     = [3.5,5.5,7.5];
            y.plot = 0;
                        
            %  > Plot variables. 
            Fig_Tools.Map_2D(fig,msh,y);
            if ~msh.f.logical(f) && inp.m.cls
                n    = ceil(max(p)./2);
                Q    = A3_2D.Q_1D_2(n);
                x_ip = xu(Q.x,msh.f.xy.v{f});
                Fig_Tools.Var_2D_2(x_ip,fig.C(1,:),"o",ms(2).*Q.w);
            end
            Fig_Tools.Var_2D_1(msh.f.xy.v  (f),fig.C(1,:),"-",lw);
            Fig_Tools.Var_2D_2(msh.c.c.xy.c   ,"k","o",ms(1));
            for i = 1:numel(sc)
                Fig_Tools.Var_2D_2(msh.c.c.xy.c(sc{i},:),fig.C(i,:),"s",ms(3));
                Fig_Tools.Var_2D_4(msh.c.c.xy.v(sc{i})  ,fig.C(i,:)    ,fig.fa);
            end
            for i = 1:numel(sf)
                Fig_Tools.Var_2D_2(msh.f.xy.c  (sf{i},:),fig.C(i,:),"s",ms(3));
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        function [] = Plot_2(fig,msh,V,y)
            %  > Auxiliary variables.
            fa = 1;
            
            %  > Plot variables.
            Fig_Tools.Var_2D_4(msh.c.c.xy.v,V,fa);
            Fig_Tools.Map_2D  (fig,msh,y);
        end
        %  > 2.2.2. -------------------------------------------------------
        function [y] = Set_y(V,flag)
            for i = 1:size(V,2)
                %  > "lim".
                switch flag
                    case 0, [y{i}.lim(1),y{i}.lim(2)] = MinMaxElem(V(:,i));
                    case 1, [y{i}.lim(1),y{i}.lim(2)] = MinMaxElem(V);
                    otherwise
                        return;
                end
                %  > "title".
                switch i
                    case 1, y{i}.title = join(["$|\bar{\tau}_{c^{\left(a\right)}}|$"]);
                    case 2, y{i}.title = join(["$|e_{c^{\left(a\right)}}|$"]);
                    otherwise
                        return;
                end
                y{i}.plot = 1;
            end
        end
    end
end