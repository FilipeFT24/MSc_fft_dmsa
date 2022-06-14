classdef Plot_2D_2
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,msh,obj)
            %  > Auxiliary variables.
            x.e = 0; %  > Export.
            x.l = 2; %  > Label.
            x.z = 0; %  > Zoom.
            fig = Fig_Tools.Set_fig(1,x);
            
            %  > Select variables.
            if inp.p.t == 2 && inp.plot{2}(1)
                n = 1;
                for i = 1:size(obj,2)
                    V  {1}(i,1) = obj(i).e.a.n_abs.t.f(1,3); %  > \tau_f.
                    V  {1}(i,2) = obj(i).e.a.n_abs.c    (1); %  >    e_c.
                    V  {2}(i,1) = obj(i).e.a.n_abs.t.f(3,3); %  > \tau_f.
                    V  {2}(i,2) = obj(i).e.a.n_abs.c    (3); %  >    e_c.
                    NNZ   (i,1) = obj(i).m{n}.nnz.At;
                end
                L1 = Plot_2D_2.Set_Y([0,1]);
                L3 = Plot_2D_2.Set_Y([0,3]);
            end
            %  > Plot variables.
            if inp.p.t == 2 && inp.plot{2}(1) 
                figure; set(gcf,'Units','pixels','Position',fig.pos);
                subplot(1,2,1);
                Plot_2D_2.Plot_1(fig,NNZ,V{1},L1);
                subplot(1,2,2);
                Plot_2D_2.Plot_1(fig,NNZ,V{2},L3);
                if x.e
                    Fig_Tools.SubPlot_pdf(fig.dir,'Plot_2D_2(1).pdf');
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 2.1.1. -------------------------------------------------------
        function [] = Plot_1(fig,NNZ,V,L)
            %  > Auxiliary variables.
            fig.fid = "Plot_2D_1(1)";
            y.dY    = [1,0];
            y.nc    = 2;
            y.plot  = 1;
            
            %  > Plot variables.
            [y.L,y.P,y.lim] = Fig_Tools.Var_1D_1(fig,repmat("-o",size(V,1),1),L,NNZ,V);         
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
                        L{2} = join(["$\|e_{c_{\left(a\right)}}\|_{_{",S1(1),"}}$"]);
                otherwise
                    return;
            end
        end
    end
end