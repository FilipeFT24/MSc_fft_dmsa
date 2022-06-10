clear, clc, close all;
%  > ----------------------------------------------------------------------
run = 2;
%  > ----------------------------------------------------------------------
switch run
    case 1
        %% > 1D.
        %  > Data.
        rng default;
        n  = 21;
        xd = linspace(-1,1,n);
        ys = xd.^2+0.25.*rand(1,n);
        z  = 0;
        %  > Fitted data.
        [xf,yf] = Fit_1D(xd,ys',z,1);
        %  > Plot...
        hold on;
        plot(xd,ys,'-or');
        plot(xf,yf,'-ob'); xlim([-1,1]); ylim([0,2]); legend('Data','Fit');
    case 2
        %% > 2D.
        %  > Data.
        rng default;
        n       = 5;
        m       = 4;
        xd      = repmat(linspace(0.3,0.7,n),m,1);
        yd      = [zeros(1,n);repmat(linspace(0.05,0.25,m-1),n,1)'];
        z       = [0.5,0];
        Xd(:,1) = reshape(xd,[],1);
        Xd(:,2) = reshape(yd,[],1);
        ys      = [-1,1,2,3,0,1,2,3,-1,1,2,3,-1,1,2,3,-1,1,2,3]';
        %  > Fitted data.
        yf      = Fit_2D(Xd,ys,z,1);
        %  > Plot...
        hold on;
        plot(1:numel(ys),ys,'or');
        plot(1:numel(yf),yf,'ob');
        %  plot(Xd(:,1),Xd(:,2),'ob'); xlim([0,1]); ylim([0,1]);
    otherwise
        return;
end
%  > ----------------------------------------------------------------------
%% > 1D.
function [Df] = Df_1D(x,xf)
    c  = [1,1,1,1];
    e  = [0,1,2,3];
    Df = c.*(x-xf)'.^e(1,:);
end
function [xf,yf] = Fit_1D(xd,ys,z,cls)
    Df = Df_1D(xd,z);
    switch cls
        case 0
            % >> ULS.
            cf = (Df'*Df)\Df'*ys;
        case 1
            % >> CLS.
            %  > Constraints.
            xc = [-0.35,0.35];
            b  = [ 0.30,0.30];
            Cf = Df_1D(xc,z);
            %  > Coefficients.
            %  > #1.
            mf = cls_c(Df'*Df,Cf,Df,b');
            cf = mf{1}*ys+mf{2};
            %  > #2.
            df = [1,0,0,0];
            C1 = df*mf{1};
            C2 = df*mf{2};
        otherwise
            return;
    end
    %  > Dd (points to evaluate polynomial at).
    switch cls
        case 0, xf = xd;
        case 1, xf = sort([xd,xc]); 
        otherwise
            return;
    end
    Dd = Df_1D(xf,z);
    yf = Dd*cf;
end
function [mf] = cls_c(DTD,C,D,b)
    %  > Auxiliary variables.
    C_DTD_CT = C*(DTD\C');
    C_DTD_DT = C*(DTD\D');
    %  > Terms.
    mf{1}    = DTD\(D'-C'*(C_DTD_CT\C_DTD_DT));
    mf{2}    = DTD\(C'*(C_DTD_CT\b));
end
%  > ----------------------------------------------------------------------
%% > 2D.
function [Df] = Df_2D(Xd,z,type)
    switch type
        case 1
            %  > \phi_f.
            c  = [1,1,1,1,1,1,1,1,1,1];
            e  = [0,1,2,3,0,1,2,0,1,0;...
                  0,0,0,0,1,1,1,2,2,3];
        case 2
            %  > \nabla\phi_f(x).
            c  = [0,1,2,3,0,1,2,0,1,0];
            e  = [0,0,1,2,0,0,1,0,0,0;...
                  0,0,0,0,1,1,1,2,2,3];
    end
    Df = c.*(Xd(:,1)-z(1)).^e(1,:).*(Xd(:,2)-z(2)).^e(2,:);
end
function [yf] = Fit_2D(Xd,ys,z,cls)
    %  > Constraint(s) (location).
    Xc = [0.45,0;0.55,0];
    for i = 1:size(Xc,1)
        ic(:,i) = round(Xd(:,1),10) == round(Xc(i,1),10)...
            & round(Xd(:,2),10) == round(Xc(i,2),10);
    end
    switch cls
        case 0
            % >> ULS.
            Df = Df_2D(Xd,z,1); Tf = (Df'*Df)\Df';
            cf = Tf*ys;
        case 1
            % >> CLS.
            for i = 1:size(ic,2)
                s{i} = find(ic(:,i));
            end
            rem = [s{:}]; Xd(rem,:) = []; ys(rem) = [];
            Df = Df_2D(Xd,z,1);
            b  = [-0.35,-0.35];
            Cf = Df_2D(Xc,z,1);
            %  > Coefficients.
            %  > #1.
            mf = cls_c(Df'*Df,Cf,Df,b');
            cf = mf{1}*ys+mf{2};
            %  > #2.
            df = Df_2D(Xc,z,2);
            C1 = sum(df*mf{1},1)./2;
            C2 = sum(df*mf{2},1)./2;
            
            ll=1;
            
        otherwise
            return;
    end
    %  > Dd (points to evaluate polynomial at).
    switch cls
        case 0, Xf = Xd;
        case 1, Xf = sort([xd,xc]); 
        otherwise
            return;
    end
    Dd = Df_2D(Xf,z);
    yf = Dd*cf;
end