classdef Fig_1_1D
    methods (Static)
        %% > Fig #1.
        function [] = Plot(msh,pde)
            %  > Auxiliary variables.
            Exp        = false;
            fig        = Fig_Tools_1D.Set_fig(Exp);
            fig.Folder = "../[Figures]/[1D]/Fig_2";
            N          = [1,2,3];

            if ~Exp
                % >> #1.
                figure(N(1)); set(gcf,'Units','pixels','Position',[150,100,1250,600]);
                %  > #1.
                subplot(1,2,1);
                Fig_1_1D.Plot_1(msh,pde,fig,Exp);
                %  > #2.
                subplot(1,2,2);
                Fig_1_1D.Plot_2(msh,pde,fig,Exp);
            else
                % >> #1.
                figure(N(1)); set(gcf,'Units','pixels','Position',[350,100,850,600]);
                Fig_1_1D.Plot_1(msh,pde,fig,N(1),Exp);
                % >> #2.
                figure(Fig(2)); set(gcf,'Units','pixels','Position',[350,100,850,600]);
            end
        end
        % >> 1. -----------------------------------------------------------
        function [] = Plot_1(msh,pde,fig,Exp)
            %  > Auxiliary variables (colors/labels).
            C    = linspecer(9,'qualitative');
            L{1} = "$|\tau_{f}^{\phantom{\nabla}\phi}|$";
            L{2} = "$|\tau_{f}^{\nabla\phi}|$";
            L{3} = "$|\tau_{f}^{\phantom{\nabla\phi}}|$";
            L{4} = "$|\!|\tau_{f}^{\phantom{\nabla}\phi}|\!|_{1}$";
            L{5} = "$|\!|\tau_{f}^{\nabla\phi}|\!|_{1}$";
            L{6} = "$|\!|\tau_{f}^{\phantom{\nabla\phi}}|\!|_{1}$";
            L{7} = "$\textrm{Error magnitude}$";
             
            %  > Plot variables.
            [P1,Y1] = Fig_Tools_1D.Var_1(fig,C,msh.f.Xv,pde.e.t.f_abs);
            [P2,Y2] = Fig_Tools_1D.Var_2(fig,C,msh.f.Xv,pde.e.t.n_abs.f);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,L,[P1,P2],[Y1;Y2]); 
            %  > Export(?).
            if Exp
                Fig_Tools_1D.Export_PDF("XXX",fig.Folder)
            end
        end
        % >> 2. -----------------------------------------------------------
        function [] = Plot_2(msh,pde,fig,Exp)
            %  > Auxiliary variables (colors/labels).
            C    = linspecer(9,'qualitative');
            L{1} = "$\phantom{|\!}|e_{c}|$";
            L{2} = "$\phantom{|\!}|\tau_{c}|$";
            L{3} = "$|\!|e_{c}|\!|_{1}$";
            L{4} = "$|\!|\tau_{c}|\!|_{1}$";
            L{5} = "$\textrm{Error magnitude}$";
             
            %  > Plot variables.
            [P1,Y1] = Fig_Tools_1D.Var_1(fig,C,msh.c.Xc,pde.e.c.c_abs);
            [P2,Y2] = Fig_Tools_1D.Var_2(fig,C,msh.c.Xc,pde.e.c.n_abs);
            [P3,Y3] = Fig_Tools_1D.Var_1(fig,C,msh.c.Xc,pde.e.t.c_abs);
            [P4,Y4] = Fig_Tools_1D.Var_2(fig,C,msh.c.Xc,pde.e.t.c_abs);
            %  > Set axis/legend,etc.
            Fig_Tools_1D.Set_Plot(fig,msh,L,[P1,P2,P3,P4],[Y1;Y2;Y3;Y4]); 
            %  > Export?
            if Exp
                Fig_Tools_1D.Export_PDF("XXX",fig.Folder)
            end
        end
    end
end