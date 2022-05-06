classdef Tools_1
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> Sort structures.
        % >> 1.1. ---------------------------------------------------------
        %  > Sort "msh" (2D) fields.
        function [msh] = Sort_msh_2D(msh)
            % >> msh.
            msh        = orderfields(msh       ,{'c','d','f','v','struct'});
            %  > c.
            msh.c      = orderfields(msh.c     ,{'c','f','h','logical','Nc','Volume'});
            msh.c.c    = orderfields(msh.c.c   ,{'nb','xy'});
            msh.c.c.nb = orderfields(msh.c.c.nb,{'f','v'});
            msh.c.c.xy = orderfields(msh.c.c.xy,{'c','v'});
            msh.c.f    = orderfields(msh.c.f   ,{'if','xy','Sf'});
            msh.c.f.xy = orderfields(msh.c.f.xy,{'c','v'});
            msh.c.h    = orderfields(msh.c.h   ,{'h','xy'});
            %  > d.
            msh.d      = orderfields(msh.d     ,{'h'});
            %  > f.
            msh.f      = orderfields(msh.f     ,{'ic','iv','logical','Nf','xy'});
            msh.f.xy   = orderfields(msh.f.xy  ,{'c','v'});
            %  > v.
            msh.v      = orderfields(msh.v     ,{'ic','if','logical'});
        end
        
        %% > 2. -----------------------------------------------------------
        % >> Other functions.
        % >> 2.1. ---------------------------------------------------------
        %  > Similar (faster) version of other built-in functions (to compute distance between 2 points).
        function [D] = dist(XY)
            D = sqrt((XY(1,1)-XY(2,1)).^2+(XY(1,2)-XY(2,2)).^2);
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Similar (faster) version of built-in function "mean".
        function [y] = mean(x,j)
            y = sum(x,j,'default','includenan')./size(x,j);
        end
        
        %% > 3. -----------------------------------------------------------
        %  > Set V=(x,y).
        % >> 3.1. ---------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        function [Xd,Yd] = msh_xy(h,t,XLim,YLim)
            %  > Limits.
            [XM(1),XM(2)] = MinMaxElem(XLim); dX = XM(2)-XM(1); Nd(1) = round(dX./h);
            [YM(1),YM(2)] = MinMaxElem(YLim); dY = YM(2)-YM(1); Nd(2) = round(dY./h);
            
            %  > Select #example...
            switch t
                case 1
                    %  > #1.
                    [Xd,Yd] = Tools_1.E1(XM,YM,Nd);
                case 2
                    %  > #2.
                    [Xd,Yd] = Tools_1.E2(XM,YM,Nd);
                case 3
                    %  > #3.
                    [Xd,Yd] = Tools_1.E3(XM,YM,Nd);
                case 4
                    %  > #4.
                    [Xd,Yd] = Tools_1.E4(XM,YM,Nd);
                case 5
                    %  > #5.
                    [Xd,Yd] = Tools_1.E5(XM,YM,Nd);
                case 6
                    %  > #5.
                    [Xd,Yd] = Tools_1.E1(XM,YM,Nd);
                otherwise
                    return;
            end  
        end
        %  > 3.1.1. -------------------------------------------------------
        %  > Uniform distribution.
        function [Xd,Yd] = E1(XM,YM,Nd)
            %  > Vd.
            Vd{1}   = linspace(XM(1),XM(2),Nd(1)+1);
            Vd{2}   = linspace(YM(1),YM(2),Nd(2)+1);
            %  > X/Yd.
            [Xd,Yd] = meshgrid(Vd{1},Vd{2});
        end
        %  > 3.1.2. -------------------------------------------------------
        %  > Transformation #1: x(u,v) = u+a*sin(2*pi*u)*sin(2*pi*v).
        %                       y(u,v) = v+a*sin(2*pi*u)*sin(2*pi*v).
        function [Xd,Yd] = E2(XM,YM,Nd)
            %  > Auxiliary variables.
            c (1)   = 2.*pi;
            c (2)   = 0.5./c(1);
            %  > Vd.
            Vd{1}   = linspace(XM(1),XM(2),Nd(1)+1);
            Vd{2}   = linspace(YM(1),YM(2),Nd(2)+1);
            %  > X/Yt.
            [Xt,Yt] = meshgrid(Vd{1},Vd{2});
            %  > X/Yd.
            f {1}   = @(u,v) u+c(2).*sin(c(1).*u).*sin(c(1).*v);
            f {2}   = @(u,v) v+c(2).*sin(c(1).*u).*sin(c(1).*v);
            Xd      = f{1}(Xt,Yt);
            Yd      = f{2}(Xt,Yt);
        end
        %  > 3.1.3. -------------------------------------------------------
        %  > Transformation #2: x(u,v) = L*(i-1)/(im-1)*[a*(im-i)/(im-1)*sin(w*H/L*(j-1)/(jm-1))+1].
        %                       y(u,v) = H*(j-1)/(jm-1)*[a*(jm-j)/(jm-1)*sin(w*H/L*(i-1)/(im-1))+1].
        function [Xd,Yd] = E3(XM,YM,Nd)
            %  > Auxiliary variables.
            L       = dX;
            H       = dY;
            a       = 0.1;
            w       = 10.*pi;
            %  > i/jm.
            [im,jm] = meshgrid(1:Nd(1)+1,1:Nd(2)+1);
            %  > X/Yd.
            f {1}   = @(i,j) L.*(i-1)./Nd(1).*(a.*(Nd(1)+1-i)./Nd(1).*sin(w.*H./L.*(j-1)./Nd(2))+1);
            f {2}   = @(i,j) H.*(j-1)./Nd(2).*(a.*(Nd(2)+1-j)./Nd(2).*sin(w.*H./L.*(i-1)./Nd(1))+1);
            Xd      = f{1}(im,jm);
            Yd      = f{2}(im,jm);
        end
        %  > 3.1.4. -------------------------------------------------------
        %  > Transformation #3 (from Artur's thesis).
        function [Xd,Yd] = E4(XM,YM,Nd)
            % > Displacement percentage (0<g<0.5).
            g = 0.25;
            %  > Check...
            if g >= 0.5
                return;
            end
            %  > Vd.
            for i = 1:numel(Nd)
                if rem(Nd(i),2)
                    h_ref(i) = 1./(Nd(i)+2.*g);
                else
                    h_ref(i) = 1./Nd(i);
                end
                for j = 1:Nd(i)+1
                    switch j
                        case 1
                            Vd{i}(j) = 0;
                        otherwise
                            if ~rem(j,2)
                                Vd{i}(j) = Vd{i}(j-1)+h_ref(i).*(1+2.*g);
                            else
                                Vd{i}(j) = Vd{i}(j-1)+h_ref(i).*(1-2.*g);
                            end
                    end
                end
            end
            %  > X/Yd.
            [Xd,Yd] = meshgrid(Vd{1},Vd{2});
        end
        %  > 3.1.5. -------------------------------------------------------
        %  > Transformation #1: x(u,v) = u+g/pi*sin(pi*u).
        %                       y(u,v) = v+g/pi*sin(pi*v).
        function [Xd,Yd] = E5(XM,YM,Nd)
            %  > Auxiliary variables.
            g       = 0.95;
            %  > Vd.
            Vd{1}   = linspace(XM(1),XM(2),Nd(1)+1);
            Vd{2}   = linspace(YM(1),YM(2),Nd(2)+1);
            %  > X/Yt.
            [Xt,Yt] = meshgrid(Vd{1},Vd{2});
            %  > X/Yd.
            f {1}   = @(u) u;
            f {2}   = @(v) v-g/pi*sin(pi*v);
            Xd      = f{1}(Xt);
            Yd      = f{2}(Yt);
        end

        %% > 4. -----------------------------------------------------------
        %  > Analytic functions.
        % >> 4.1. ---------------------------------------------------------
        function [fh] = func(c,f_type)
            %  > Symbolic variables.
            syms x y;
            
            switch f_type
                case 1
                    xc = 0.5;
                    yc = 0.5;
                    i  = 100;
                    f  = exp(-i.*((x-xc).^2+(y-yc).^2));
                otherwise
                    return;
            end
            fh = matlabFunction(f,'Vars',{x,y});
        end
    end
end