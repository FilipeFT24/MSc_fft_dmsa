classdef Plot_2D_2
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,msh,obj,v)
            %  > Auxiliary variables.
            x  (1).e = 0;
            x  (1).l = 2;
            x  (1).z = 0;
            fig{1}   = Fig_Tools.Set_fig(1,x(1));
            x  (2).e = 0;
            x  (2).l = 3;
            x  (2).z = 0;
            fig{2}   = Fig_Tools.Set_fig(0,x(2));
            
            %  > Select variables.
            if inp.T == 2 && inp.Plot{2}(1)
                for i = 1:size(obj,2)
                    Y{1}(i,1) = obj(i).e.a.n_abs.t.f(1,3); %  > \tau_f.
                    Y{1}(i,2) = obj(i).e.a.n_abs.t.c(1);   %  > \tau_c.
                    Y{1}(i,3) = obj(i).e.a.n_abs.c  (1);   %  >     ec.
                    Y{2}(i,1) = obj(i).e.a.n_abs.t.f(3,3); %  > \tau_f.
                    Y{2}(i,2) = obj(i).e.a.n_abs.t.c(3);   %  > \tau_c.
                    Y{2}(i,3) = obj(i).e.a.n_abs.c  (3);   %  >     ec.
                    X   (i,1) = obj(i).m.nnz.At;
                end
                L{1} = Plot_2D_2.Set_Y([0,1]);
                L{2} = Plot_2D_2.Set_Y([0,3]);
            end
            if inp.T == 2 && inp.Plot{2}(2), m = size(obj,2); n = numel(obj(m).s);
                for i = 1:n
                    e(:,i)     = obj(m).e.a.t.f_abs(:,end);
                    p(:,i)     = obj(m-1).s.u.p      (:,i);
                    s  (i).r.f = v  (m-1).r.f;
                    s  (i).r.s = v  (m-1).r.s;
                end
                nnz_A = obj(m).m.nnz.A;
            end
            %  > Plot variables.
            if inp.T == 2 && inp.Plot{2}(1)
                fprintf("Plotting (2/1)...\n");
                figure;
                for i = 1:2
                    subplot(1,2,i); Plot_2D_2.Plot_1(fig{1},X,Y{i},L{i});
                end
                if x(1).e
                    exportgraphics(gcf,'Plot_2D_2(1).pdf','ContentType','Vector');
                end
            end
            if inp.T == 2 && inp.Plot{2}(2)
                fprintf("Plotting (2/2)...\n");
                figure;
                for i = 1:n
                    subplot(1,n,i); Plot_2D_2.Plot_2(fig{2},msh,nnz_A,e(:,i),p(:,i),s(i));
                end
                if x(2).e
                    exportgraphics(gcf,'Plot_2D_2(2).pdf','ContentType','Vector');
                end
            end
            %  > Move to...
            if inp.Plot{2}(1) && x(1).e
                movefile('Plot_2D_2(1).pdf','[Post-processing]/[.PDF Files]');
            end
            if inp.Plot{2}(2) && x(2).e
                movefile('Plot_2D_2(2).pdf','[Post-processing]/[.PDF Files]');
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 2.1.1. -------------------------------------------------------
        function [] = Plot_1(fig,X,Y,L)
            %  > Auxiliary variables.
            Z.DY   = [1,0];
            Z.NC   = 1;
            Z.Plot = 1;
            
            %  > Plot variables.
            [Z.L,Z.P,Z.lim] = Fig_Tools.Var_1D_1(fig,["-v","-s","-o"],L,X,Y);         
            if Z.lim (1) < fig.TrshV
                Z.lim(1) = fig.TrshV;
            end
            Fig_Tools.Map_1D(fig,X,Z);        
        end
        %  > 2.1.2. -------------------------------------------------------
        function [L] = Set_Y(opt)
            S(1) = Fig_Tools.S2(opt(2));
            switch opt(1)
                %  > w/o predicted.
                case 0, L{1} = join(["$\|\bar{\tau}_{f_{\left(a\right)}}^{\phi}\|_{_{",S(1),"}}$"]);
                        L{2} = join(["$\|\bar{\tau}_{c_{\left(a\right)}}\|_{_{",S(1),"}}$"]);
                        L{3} = join(["$\|e_{c_{\left(a\right)}}\|_{_{",S(1),"}}$"]);
                otherwise
                    return;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        function [] = Plot_2(fig,msh,nnz,e,p,v)
            %  > Auxiliary variables.
            for i = 1:msh.f.Nf
                Length(i,1) = sqrt(diff(msh.f.xy.v{i}(:,1)).^2+diff(msh.f.xy.v{i}(:,2)).^2);
            end
            MS = Fig_Tools.Convert_MS(mean(Length).*1.05);
            T  = 0;

            %  > Plot variables.
            Fig_Tools.Map_2D(fig,msh);
            for i = 1:msh.f.Nf
                Fig_Tools.Var_2D_2(msh.f.xy.c(i,:),fig.C(ceil(p(i)./2),:),"o",MS);
            end
            if T
                for i = 1:msh.f.Nf
                    text(msh.f.xy.c(i,1),msh.f.xy.c(i,2),cellstr(num2str(i)),'Fontsize',fig.FS{2}./2);
                end
            end
            Fig_Tools.Var_2D_3(msh.f.xy.c(v.r.f,:),"k",1.00,"o",1.15.*MS); %  > "f" (flagged).
            Fig_Tools.Var_2D_3(msh.f.xy.c(v.r.s,:),"k",1.00,"s",0.85.*MS); %  > "s" (selected).
            Fig_Tools.Var_2D_2(msh.f.xy.c(e >= max(e)-10.^(ceil(log10(max(e))-1)-3),:),"k","+",0.25.*MS);
            %  > Legend.
            Plot_2D_2.Set_L(fig,MS,p);
        end
        %  > 2.2.1.1. -----------------------------------------------------
        function [] = Set_L(fig,MS,p)
            j = ceil(RunLength(sort(p))./2);
            for i = 1:numel(j)
                P{i} = plot(NaN,NaN,"o",'Color',fig.C(j(i),:),'MarkerFaceColor',fig.C(j(i),:),'MarkerSize',MS);
                L{i} = num2str(j(i));
            end
            A = get(gca,'Position');
            L = legend([P{:}],[L],'Location','NortheastOutside','FontSize',fig.FS{3},'NumColumns',1);
            set  (get(L,'title'),'String','Level','FontSize',fig.FS{3});
            set  (gca,'Position',A);
        end
    end
end