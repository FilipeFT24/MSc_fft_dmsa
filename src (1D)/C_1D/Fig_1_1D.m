classdef Fig_1_1D
    methods (Static)         
        %% > Wrap-up Fig_1 (1D).
        function [] = WrapUp_Fig_1_1D(Plot_1,Exp_1,Fig,msh)
            if Plot_1
                %  > Select...
                non_empty = find(cellfun(@isempty,msh.s.f));
                iF        = non_empty(randperm(length(non_empty),1));               
                %  > Figure.
                figure(Fig); set(gcf,'Units','pixels','Position',[250,150,1000,500]);
                Fig_1_1D.Plot(msh,msh.f.NF);
                %  > Export as .pdf.
                if Exp_1
                    Fig_Tools_1D.Export_PDF('Fig_1','../[Figures]/[1D]/Fig_1');
                end
            end
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [] = Plot(msh,iF)
            C  = linspecer(9,'qualitative');
            X1 = msh.s.xt{iF};
            X2 = setdiff(msh.f.Xv,X1);
            X3 = [setdiff(msh.f.Xv,msh.f.Xv(iF));repelem(0,msh.f.NF-1)];
            
            hold on;
            plot(X3(1,:),X3(2,:),'-k','MarkerSize',10.0);
            plot(X3(1,:),X3(2,:),'|','Color','k','MarkerFaceColor','k','MarkerSize',10.0);
            plot(X1,0,'s','Color',C(2,:),'MarkerFaceColor',C(2,:),'MarkerSize',5.0);
            plot(msh.f.Xv(iF),0,'|','Color',C(1,:),'MarkerFaceColor',C(1,:),'MarkerSize',10.0,'Linewidth',1.5);
            %  > Axis.
            Fig_Tools_1D.ChangeLook_1D(false,false,msh.f.Xv,10,"$x$","$y$",20,12);
            ax = gca; 
            ax.XAxis.Visible = 'off'; ax.YAxis.Visible = 'off';
        end
    end
end      