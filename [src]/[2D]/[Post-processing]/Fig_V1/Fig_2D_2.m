classdef Fig_2D_2
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,msh,obj_e)
            %  > Auxiliary variables.
            exp   = 0;
            run   = 0;
            zoom  = 0;
            fig   = Fig_Tools.Set_fig(exp,run,zoom);
            n     = 1; %  > {n}. 
            inp_m = A1_2D.Set_msh(msh.d.h);
            flag  = inp_m.p == "s" && inp_m.t == 0;
        
            if ~exp
                if inp.plot(2)
                    %  > Colormap.
                    figure; set(gcf,'Units','pixels','Position',fig.Position);
                    subplot(1,2,1);
                    Fig_2D_2.Plot_1(msh,1,obj_e.a.t.c_abs,fig);
                    subplot(1,2,2);
                    Fig_2D_2.Plot_1(msh,2,obj_e.a.c.c_abs,fig);
                    %  > 1D.
                    if flag
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         figure; set(gcf,'Units','pixels','Position',fig.Position);
%                         j = 1; %  > (:,j).
%                         subplot(1,2,j);
%                         Fig_2D_2.Plot_2(inp_m,msh,j,obj.e.a{n}.t.c_abs,obj.e.a{n}.c.c_abs,fig);
%                         j = 2; %  > (:,j).
%                         subplot(1,2,j);
%                         Fig_2D_2.Plot_2(inp_m,msh,j,obj.e.a{n}.t.c_abs,obj.e.a{n}.c.c_abs,fig);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         figure; set(gcf,'Units','pixels','Position',fig.Position);
%                         j = 1; %  > (:,j).
%                         subplot(1,2,j);
%                         Fig_2D_2.Plot_3(inp_m,msh,j,obj.e.a{n}.t.f_abs.y,fig);
                    end
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [] = Plot_1(msh,s,v,fig)
            %  > Auxiliary variables.
            fig.fid         = "2D_2";
            fig.LW          = 0.1;
            fig.FT_1        = 12.00;
            y.str           = 'thermal';
            y.title         = Fig_2D_2.Set_Legend_1(s);
            [y.c(1),y.c(2)] = MinMaxElem(v);
            
            %  > Plot variables.
            Fig_Tools.Var_2D_3(msh.c.c.xy.v,v,1);            
            %  > Axis/legend,etc.
            Fig_Tools.Map_2D_1(msh,fig);
            Fig_Tools.Map_2D_2(fig,y);
            Fig_Tools.Map_2D_3(fig);
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_2(inp_m,msh,j,V1,V2,fig)
            %  > Find faces...
            d = Tools_1.setdiff(1:size(msh.c.c.xy.c,2),j);
            c = RunLength(sort(msh.c.c.xy.c(:,d)));
            k = msh.d.h./10E2;
            
            %  > Select variables.
            for i = 1:numel(c)
                l (:,i) = find(msh.c.c.xy.c(:,d) > c(i)-k & msh.c.c.xy.c(:,d) < c(i)+k);
                xy(:,i) = msh.c.c.xy.c(l(:,i),j);
                U1(:,i) = V1(l(:,i));
                U2(:,i) = V2(l(:,i));
            end
            %  > Plot variables.
            L1{1} = Fig_2D_2.Set_Legend_1(1);
            L1{2} = Fig_2D_2.Set_Legend_1(2);
            M1    = ["--o",":v"];
            for i = 1:numel(c)
                [L1{i},P1,YX] = Fig_Tools.Var_1D_1(fig,M1,L1{1},xy(:,i),[U1(:,i)]);
                if i == 1
                    Y1 = YX;
                else
                    Y1 = min(YX,Y1);
                end
            end
            %  > Axis/legend,etc.
            Fig_Tools.Map_1D_2(fig,L1,P1,0,inp_m.Lim(j,:),Y1,[-1,1],2); Fig_2D_2.Set_Legend_2(fig,j);
        end
        % >> 2.2. ---------------------------------------------------------
        function [] = Plot_3(inp_m,msh,j,V1,fig)
            %  > Find faces...
            d = Tools_1.setdiff(1:size(msh.f.xy.c,2),j);
            c = 0:0.2:1;%RunLength(sort(msh.f.xy.c(:,d)));
            k = msh.d.h./10E2;
            
            %  > Select variables.
            for i = 1:numel(c)
                l (:,i) = find(msh.f.xy.c(:,d) > c(i)-k & msh.f.xy.c(:,d) < c(i)+k);
                xy(:,i) = msh.f.xy.c(l(:,i),j);
                U1(:,i) = V1(l(:,i),3);
            end
            hold on;
            for i = 1:numel(c)
                plot(xy(:,i),U1(:,i));
            end
            
            
            
        end
        % >> 2.4 ----------------------------------------------------------
        %  > 2.4.1 --------------------------------------------------------
        function [L] = Set_Legend_1(t)
            switch t
                case 1, L = join(["$|\bar{\tau}_{c^{\left(a\right)}}|$"]);
                case 2, L = join(["$|e_{c^{\left(a\right)}}|$"]);
                otherwise
                    return;
            end
        end
        %  > 2.4.2 --------------------------------------------------------
        function [] = Set_Legend_2(fig,j)
            switch j
                case 1
                    L{1} = "$y$";
                    L{2} = "$\textrm{Error magnitude}$";
                case 2
                    L{1} = "$x$";
                    L{2} = "$\textrm{Error magnitude}$";
            end
            xlabel(L{1},'FontSize',fig.FT_3(1),'Interpreter','latex');
            ylabel(L{2},'FontSize',fig.FT_3(1),'Interpreter','latex');
        end
    end
end