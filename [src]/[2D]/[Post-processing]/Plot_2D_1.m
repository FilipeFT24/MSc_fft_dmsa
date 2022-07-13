classdef Plot_2D_1
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,msh,obj)
            %  > Auxiliary variables.
            x.e = 1;
            x.l = 3;
            x.z = 0;
            fig = Fig_Tools.Set_fig(0,x);
            
            %  > Select variables.
            if inp.Plot{1}(1)
                f = [1,2];
            end
            if inp.Plot{1}(2)
                V = [obj(end).e.a.t.c_abs,obj(end).e.a.c.c_abs];
            end
            %  > Plot variables.
            if inp.Plot{1}(1)
                fprintf("Plotting (1/1)...\n");
                figure;
                for i = 1:2
                    subplot(1,2,i); Plot_2D_1.Plot_1(fig,inp,msh,f(i),obj(end).s.sc{f(i)},obj(end).s.sf{f(i)},obj(end).s.u.p(f(i),:),obj(end).f.qd.xu);
                end
                if x.e
                    exportgraphics(gcf,'Plot_2D_1(1).pdf','ContentType','Vector');
                end
            end
            if inp.Plot{1}(2)
                fprintf("Plotting (1/2)...\n");
                figure;
                Z = Plot_2D_1.Set_Z(fig,1,V);
                for i = 1:2
                    subplot(1,2,i); Plot_2D_1.Plot_2(fig,msh,x.e,V(:,i),Z{i});
                end
                if x.e
                    exportgraphics(gcf,'Plot_2D_1(2).pdf','ContentType','Vector');
                end
            end
            %  > Move to...
            if inp.Plot{1}(1) && x.e
                movefile('Plot_2D_1(1).pdf','[Post-processing]/[.PDF Files]');
            end
            if inp.Plot{1}(2) && x.e
                movefile('Plot_2D_1(2).pdf','[Post-processing]/[.PDF Files]');
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1(fig,inp,msh,f,sc,sf,p,xu)
            %  > Auxiliary variables.
            [x,w] = Plot_2D_1.x_ip(inp,msh,f,p,xu);
            for i = 1:msh.f.Nf
                Length(i,1) = sqrt(diff(msh.f.xy.v{i}(:,1)).^2+diff(msh.f.xy.v{i}(:,2)).^2);
            end
            MS = Fig_Tools.Convert_MS(mean(Length)./2);
            
            %  > Plot variables.
            Fig_Tools.Map_2D  (fig,msh);
            Fig_Tools.Var_2D_6(msh.c.c.xy.c,msh.c.c.xy.v,msh.f.xy.v(f),msh.f.xy.c,sc,sf,x,w,fig.C,fig.FA,fig.LW,["o","s"],[MS./2,MS]);
        end
        % >> 2.1.1. -------------------------------------------------------
        function [x,w] = x_ip(inp,msh,f,p,xu)
            if ~msh.f.logical(f) && inp.M.Cls
                Q = A3_2D.Q_1D_2(ceil(max(p)./2));
                x = xu(Q.x,msh.f.xy.v{f});
                w = Q.w./2;
            else
                x = msh.f.xy.c(f,:);
                w = 1;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        function [] = Plot_2(fig,msh,e,V,Z)
            Z.e = e;
            Fig_Tools.Var_2D_4(msh.c.c.xy.v,V,1.00);
            Fig_Tools.Map_2D  (fig,msh,Z);
        end
        %  > 2.2.2. -------------------------------------------------------
        function [Z] = Set_Z(fig,opt,V)
            for i = 1:size(V,2)
                %  > "Lim".
                switch opt
                    case 0, [Z{i}.Lim(1),Z{i}.Lim(2)] = MinMaxElem(V(:,i));
                    case 1, [Z{i}.Lim(1),Z{i}.Lim(2)] = MinMaxElem(V);
                    otherwise
                        return;
                end
                if  Z{i}.Lim(1) >= fig.TrshV
                    Z{i}.Lim(1)  = 10.^(ceil(log10(Z{i}.Lim(1)))-1);
                    Z{i}.Lim(2)  = 10.^(ceil(log10(Z{i}.Lim(2)))+0);
                else
                    Z{i}.Lim(1)  = fig.TrshV;
                    Z{i}.Lim(2)  = 10.^(ceil(log10(Z{i}.Lim(2)))+0);
                end
                %  > "Title".
                switch i
                    case 1, Z{i}.Title = join(["$|\bar{\tau}_{c^{\left(a\right)}}|$"]);
                    case 2, Z{i}.Title = join(["$|e_{c^{\left(a\right)}}|$"]);
                    otherwise
                        return;
                end
            end
        end
    end
end