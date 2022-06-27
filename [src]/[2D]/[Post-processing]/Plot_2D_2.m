classdef Plot_2D_2
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,msh,obj)
            %  > Auxiliary variables.
            x  (1).e = 0; %  > Export.
            x  (1).l = 2; %  > Label.
            x  (1).z = 0; %  > Zoom.
            fig{1}   = Fig_Tools.Set_fig(1,x(1));
            x  (2).e = 0; %  > Export.
            x  (2).l = 3; %  > Label.
            x  (2).z = 0; %  > Zoom.
            fig{2}   = Fig_Tools.Set_fig(0,x(2));
            
            %  > Select variables.
            if inp.p.t == 2 && inp.plot{2}(1)
                n = 1;
                for i = 1:size(obj,2)
                    V  {1}(i,1) = obj(i).e.a.n_abs.t.f(1); %  > \tau_f.
                    V  {1}(i,2) = obj(i).e.a.n_abs.t.c(1); %  > \tau_c.
                    V  {1}(i,3) = obj(i).e.a.n_abs.c  (1); %  >    e_c.
                    V  {2}(i,1) = obj(i).e.a.n_abs.t.f(3); %  > \tau_f.
                    V  {2}(i,2) = obj(i).e.a.n_abs.t.c(3); %  > \tau_c.
                    V  {2}(i,3) = obj(i).e.a.n_abs.c  (3); %  >    e_c.
                    NNZ   (i,1) = obj(i).m{n}.nnz.A;
                end
                L1 = Plot_2D_2.Set_Y([0,1]);
                L3 = Plot_2D_2.Set_Y([0,3]);
            end
            if inp.p.t == 2 && inp.plot{2}(2)
                e{1} = obj(end-1).e.a.t.f_abs;
                p{1} = obj(end-1).p;
                s{1} = obj(end-1).s;
                e{2} = obj(end)  .e.a.t.f_abs;
                p{2} = obj(end)  .p;
                s{2} = obj(end)  .s;
            end
            %  > Plot variables.
            if inp.p.t == 2 && inp.plot{2}(1) 
                figure; set(gcf,'Units','pixels','Position',fig{1}.pos);
                subplot(1,2,1);
                Plot_2D_2.Plot_1(fig{1},NNZ,V{1},L1);
                if x(1).z
                    zp = BaseZoom(); zp.plot;
                end
                subplot(1,2,2);
                Plot_2D_2.Plot_1(fig{1},NNZ,V{2},L3);
                if x(1).z
                    zp = BaseZoom(); zp.plot;
                end
                if x(1).e
                    Fig_Tools.SubPlot_pdf(fig{1}.dir,'Plot_2D_2(1).pdf');
                end
            end
            %  > Plot variables.
            if inp.p.t == 2 && inp.plot{2}(2) 
                figure; set(gcf,'Units','pixels','Position',fig{1}.pos);
                subplot(1,2,1);
                Plot_2D_2.Plot_2(fig{2},inp,msh,e{1},p{1},s{1},[0,size(obj,2)-2]);
                if x(2).z
                    zp = BaseZoom();
                    zp.plot;
                end
                subplot(1,2,2);
                Plot_2D_2.Plot_2(fig{2},inp,msh,e{2},p{2},s{2},[1,size(obj,2)-1]);
                if x(2).z
                    zp = BaseZoom();
                    zp.plot;
                end
                if x(2).e
                    Fig_Tools.SubPlot_pdf(fig{2}.dir,'Plot_2D_2(2).pdf');
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 2.1.1. -------------------------------------------------------
        function [] = Plot_1(fig,NNZ,V,L)
            %  > Auxiliary variables.
            y.dY   = [1,0];
            y.nc   = 1;
            y.plot = 1;
            
            %  > Plot variables.
            [y.L,y.P,y.lim] = Fig_Tools.Var_1D_1(fig,repmat("-o",size(V,2),1),L,NNZ,V);         
            if y.lim (1) < fig.trsh.v
                y.lim(1) = fig.trsh.v;
            end
            Fig_Tools.Map_1D(fig,NNZ,y);            
        end
        %  > 2.1.2. -------------------------------------------------------
        function [L] = Set_Y(flag)
            S1(1) = Fig_Tools.Set_str_2(flag(2));
            switch flag(1)
                %  > w/o predicted.
                case 0, L{1} = join(["$\|\bar{\tau}_{f_{\left(a\right)}}^{\phi}\|_{_{",S1(1),"}}$"]);
                        L{2} = join(["$\|\bar{\tau}_{c_{\left(a\right)}}\|_{_{",S1(1),"}}$"]);
                        L{3} = join(["$\|e_{c_{\left(a\right)}}\|_{_{",S1(1),"}}$"]);
                otherwise
                    return;
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        function [] = Plot_2(fig,inp,msh,e,p,s,flag)
            %  > Auxiliary variables.
            add_text = 0;
            add_cr   = 1;
            Qp       = A3_2D.Q_1D_1;
            y.plot   = 0;
            cu       = get(gca,'Units');    set(gca,'Units','Points');
            ax_p     = get(gca,'Position'); set(gca,'Units',cu);

            %  > Plot variables.
            Fig_Tools.Map_2D(fig,msh,y);
            for i = 1:msh.f.Nf
                ms_i(i,1) = sqrt(diff(msh.f.xy.v{i}(:,1)).^2+diff(msh.f.xy.v{i}(:,2)).^2)./(1.50.*diff(xlim))*ax_p(3);
            end
            ms(1) = mean(ms_i(:,1));
            ms(2) = 1.75.*ms(1);
            ms(3) = 1.50.*ms(1);
            
            for i = 1:msh.f.Nf
                if ~inp.p.iso
                    xy_p{i,1} = Qp.xu([-1/3;1/3],msh.f.xy.v{i});
                    j         = 1:size(xy_p{i},1);
                else
                    xy_p{i,1} = Qp.xu(0,msh.f.xy.v{i});
                    j         = 1;
                end
                %  > #1.
                for k = j
                    Fig_Tools.Var_2D_2(xy_p{i}(k,:),fig.C(ceil(p(i,k)./2),:),"o",ms(1));
                    if add_text
                        text(xy_p{i}(k,1),xy_p{i}(k,2),cellstr(num2str(i)),'Fontsize',5);
                    end
                end
                %  > #2.
                if add_cr
                    if ~isempty(s.c)
                        if ismembc(i,sort(s.c))
                            for k = j
                                Fig_Tools.Var_2D_3(xy_p{i}(k,:),'k',fig.lw,"s",ms(2));
                            end
                        end
                    end
                    if ~isempty(s.r)
                        if ismembc(i,sort(s.r))
                            for k = j
                                Fig_Tools.Var_2D_3(xy_p{i}(k,:),'k',fig.lw,"o",ms(3));
                            end
                        end
                    end
                end
            end
            if flag(1)
                Plot_2D_2.Set_L(fig,inp,p,"o",ms(1));
            end
            title("i = "+flag(2),'Fontsize',fig.fs{5},'Interpreter','latex');
        end
        %  > 2.2.2. -------------------------------------------------------
        function [] = Set_L(fig,inp,p,mf,ms)
            if inp.p.iso
                j = 1;
                k = ceil(RunLength(sort(p(:,j)))./2);
                for l = 1:numel(k)
                    P{l} = plot(NaN,NaN,mf,'Color',fig.C(k(l),:),'MarkerFaceColor',fig.C(k(l),:),'MarkerSize',ms);
                    L{l} = num2str(k(l));
                end
            end
            ax_kepp = get(gca,'Position');
            legend([P{:}],[L],...
                'Interpreter','latex','Location','NortheastOutside','FontSize',fig.fs{3},'NumColumns',1);
            set(gca,'Position',ax_kepp)
        end
    end
end