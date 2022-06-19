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
                    V  {1}(i,1) = obj(i).e.a.n_abs.t.f(1,3); %  > \tau_f.
                    V  {1}(i,2) = obj(i).e.a.n_abs.t.c  (1); %  > \tau_c.
                    V  {1}(i,3) = obj(i).e.a.n_abs.c    (1); %  >    e_c.
                    V  {2}(i,1) = obj(i).e.a.n_abs.t.f(3,3); %  > \tau_f.
                    V  {2}(i,2) = obj(i).e.a.n_abs.t.c  (3); %  > \tau_c.
                    V  {2}(i,3) = obj(i).e.a.n_abs.c    (3); %  >    e_c.
                    NNZ   (i,1) = obj(i).m{n}.nnz.At;
                end
                L1 = Plot_2D_2.Set_Y([0,1]);
                L3 = Plot_2D_2.Set_Y([0,3]);
            end
            if inp.p.t == 2 && inp.plot{2}(2)
                %  > Diffusive term (x) (for example).
                %    NOTE: for the conducted tests, every term uses the same stencil.
                j(1) = 2; %  > (:,j).
                k(1) = 1; %      {k}.
                e{1} = obj(end).e.a.t.f_abs(:,end);
                p{1} = obj(end).p{j(1)}{k(1)};
            end
            %  > Plot variables.
            if inp.p.t == 2 && inp.plot{2}(1) 
                figure; set(gcf,'Units','pixels','Position',fig{1}.pos);
                subplot(1,2,1);
                Plot_2D_2.Plot_1(fig{1},NNZ,V{1},L1);
                subplot(1,2,2);
                Plot_2D_2.Plot_1(fig{1},NNZ,V{2},L3);
                if x(1).e
                    Fig_Tools.SubPlot_pdf(fig{1}.dir,'Plot_2D_2(1).pdf');
                end
            end
            %  > Plot variables.
            if inp.p.t == 2 && inp.plot{2}(2) 
                figure; set(gcf,'Units','pixels','Position',fig{2}.pos);
                Plot_2D_2.Plot_2(fig{2},inp,msh,e{1},p{1});
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
            y.dY    = [1,0];
            y.nc    = 1;
            y.plot  = 1;
            
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
        function [] = Plot_2(fig,inp,msh,e,p)
            %  > Auxiliary variables.
            add_text = false;
            mf       = ["s","o"];
            ms       = 9.5;
            Qp       = A3_2D.Q_1D_1;
            y.plot   = 0;
                        
            %  > Plot variables.
            Fig_Tools.Map_2D(fig,msh,y);
            for i = 1:msh.f.Nf
                if ~inp.p.iso
                    xy_p{i,1} = Qp.xu([-1/3;1/3],msh.f.xy.v{i});
                    j         = 1:size(xy_p{i},1);
                else
                    xy_p{i,1} = Qp.xu(0,msh.f.xy.v{i});
                    j         = 1;
                end
                for k = j
                    Fig_Tools.Var_2D_2(xy_p{i}(k,:),fig.C(ceil(p(i,k)./2),:),mf(k),ms);
                    if add_text
                        %  > Auxiliary variables.
                        c     = cellstr(num2str(i));
                        shift = [-0.015,0.015];
                        %  > Plot text...
                        text(xy_p{i}(k,1)+shift(1),xy_p{i}(k,2)+shift(2),c,'Fontsize',5);
                    end
                end
            end
            Plot_2D_2.Set_L(fig,inp,p,mf,ms);
        end
        %  > 2.2.2. -------------------------------------------------------
        function [] = Set_L(fig,inp,p,mf,ms)
            if inp.p.iso
                j = 1;
                k = ceil(RunLength(sort(p(:,j)))./2);
                for l = 1:numel(k)
                    P{l} = plot(NaN,NaN,mf(j),'Color',fig.C(k(l),:),'MarkerFaceColor',fig.C(k(l),:),'MarkerSize',ms);
                    L{l} = num2str(k(l));
                end
            end
            legend([P{:}],[L],...
                'Interpreter','latex','Location','NortheastOutside','FontSize',fig.fs{3},'NumColumns',1);
        end
    end
end