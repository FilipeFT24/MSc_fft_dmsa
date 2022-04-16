classdef Fig_V1_3_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        function [] = Plot(inp,obj,msh,plot)
            %  > Auxiliary variables.
            exp = 0;
            j   = 1;
            n   = 1;
            
            if ~exp
                if ~inp.pa.adapt
                    % >> 'p-standard' run.
                    fig  = Fig_Tools_1D.Set_fig(0,exp);
                    plot = [0,1];
                    %  > #1.
                    if plot(1)
                        figure; set(gcf,'Units','pixels','Position',fig.Position);
                        subplot(1,2,1);
                        Fig_V1_1_1D.Plot_1_1(obj,msh,fig,j,n,0);
                        subplot(1,2,2);
                        Fig_V1_1_1D.Plot_1_2(obj,msh,fig,j,n,0);
                    end
                    %  > #2.
                    if plot(2)
                        figure; set(gcf,'Units','pixels','Position',fig.Position);
                        subplot(1,2,1);
                        Fig_V1_1_1D.Plot_2_1(obj.e,msh,fig,j,n,0);
                        subplot(1,2,2);
                        Fig_V1_1_1D.Plot_2_2(obj,msh,fig,j,n,0);
                    end
                else
                    % >> 'p-adaptative' run.
                    plot = [1,0];
                    %  > #1.
                    if plot(1)
                        fig = Fig_Tools_1D.Set_fig(1,exp);
                        figure; set(gcf,'Units','pixels','Position',fig.Position);
                        subplot(1,2,1);
                        Fig_V1_1_1D.Plot_3_1(inp,obj,msh,fig,n,0);
                        fig = Fig_Tools_1D.Set_fig(0,exp);
                        subplot(1,2,2);
                        Fig_V1_1_1D.Plot_3_2(inp,obj,msh,fig,n,0);
                    end
                end
            else
            end
        end

        
        % >> 3.1. ---------------------------------------------------------
        function [] = Plot_3_1(inp,obj,msh,fig,n,zoom)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');
            %fld   = ["a","p"];
            fld   = "p";
            
            %  > Set variables to plot...
            for i = 1:length(fld)
                j = fld(i);
                %  > nnz.
                for k = 1:size(obj.(j).m,1)
                    nnz.(j)(k,1) = obj.(j).m(k).nnz.Ac(1);
                    nnz.(j)(k,2) = obj.(j).m(k).nnz.Ac(2);
                    nnz.(j)(k,3) = obj.(j).m(k).nnz.At;
                end
                [x.(j)(1),x.(j)(2)] = MinMaxElem(nnz.(j));
                %  > e.
                l = 1;
                for k = 1:size(obj.(j).m,1)
                    V.(j)(k,:)   = obj.(j).e(k).(j){n}.t.n_abs.f(l,:);
                end 
            end
            [x_l(1),x_l(2)] = MinMaxElem(x.(fld));
            
      
%             M1         = ["-o","-o","-o"];
%             L1{1}      = "$\|\bar{\tau}_{f^{\left(a\right)}}^{\phantom{\nabla}\phi}\|_{1}$";
%             L1{2}      = "$\|\bar{\tau}_{f^{\left(a\right)}}^{\nabla\phi}\|_{1}$";
%             L1{3}      = "$\|\bar{\tau}_{f^{\left(a\right)}}^{\phantom{\nabla\phi}}\|_{1}$";
%             [L1,P1,Y1] = Fig_Tools_1D.Var_1(fig,M1,L1,nnz.a,V.a);
            M2         = [":v",":v",":v"];
            L2{1}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla}\phi}\|_{1}$";
            L2{2}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\nabla\phi}\|_{1}$";
            L2{3}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla\phi}}\|_{1}$";
            [L2,P2,Y2] = Fig_Tools_1D.Var_1(fig,M2,L2,nnz.p,V.p);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,L2,P2,x_l,Y2,[-1,1],2);
            %Fig_Tools_1D.Set_Plot_2(fig,[L1,L2],[P1,P2],x_l,[Y1;Y2],[-1,1],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,"X");
        end
        
        %%%%%%%%%%%%%%%%%%%
        function [] = Plot_3_2(inp,obj,msh,fig,n,zoom)
            %  > Auxiliary variables.
            fig.C = linspecer(9,'qualitative');

            V = obj.p.e(size(obj.p.e,1)).p{n}.t.f_abs(:,1:2);
            V = [V,obj.p.e(size(obj.p.e,1)).a{n}.t.f_abs(:,1:2)];
            
   
            M2         = [":v",":v","-o","-o"];
            
            L2{1}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla}\phi}\|_{1}$";
            L2{2}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\nabla\phi}\|_{1}$";
            %L2{3}      = "$\|\bar{\tau}_{f^{\left(p\right)}}^{\phantom{\nabla\phi}}\|_{1}$";
            L2{3}      = "$\|\bar{\tau}_{f^{\left(a\right)}}^{\phantom{\nabla}\phi}\|_{1}$";
            L2{4}      = "$\|\bar{\tau}_{f^{\left(a\right)}}^{\nabla\phi}\|_{1}$";
            %L2{6}      = "$\|\bar{\tau}_{f^{\left(a\right)}}^{\phantom{\nabla\phi}}\|_{1}$";
            [L2,P2,Y2] = Fig_Tools_1D.Var_1(fig,M2,L2,msh.f.Xv,V);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot_2(fig,L2,P2,msh.f.Xv,Y2,[-1,1],2);
            %Fig_Tools_1D.Set_Plot_2(fig,[L1,L2],[P1,P2],x_l,[Y1;Y2],[-1,1],2);
            %  > Export(?)/Zoom(?).
            Fig_Tools_1D.Exp_Or_Zoom(fig,zoom,"X");
        end
    end
end