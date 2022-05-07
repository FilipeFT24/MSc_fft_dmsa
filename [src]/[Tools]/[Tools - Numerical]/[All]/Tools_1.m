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
        %  > Set "struct" structure.
        % >> 3.1. ---------------------------------------------------------
        function [struct] = Set_struct(m)
            %  > ...for x/y-directions.
            for i = 1:size(m.Lim,1)
                %  > Limits.
                [L(i,1),L(i,2)] = MinMaxElem(m.Lim(i,:)); Nv(i) = round(diff(L(i,:))./m.h)+1; Nc = Nv-1;
                %  > X/Yv.
                XY_v{i} = linspace(L(i,1),L(i,2),Nv(i));
            end
            %  > X/Yt,d.
            [Xt,Yt] = meshgrid  (XY_v{1},XY_v{2});
            [Xd,Yd] = Tools_1.T1(Xt,Yt);
            
            %  > Select cell polyhedral (type)...
            switch m.p
                case "s"
                    % >> w/ squares.
                    %  > Connectivity list.
                    for i = 1:Nc(1)
                        for j = 1:Nc(2)
                            struct.ConnectivityList(Nc(2)*(i-1)+j,1) = Nv(2)*(i-1)+j;   % > SW.
                            struct.ConnectivityList(Nc(2)*(i-1)+j,2) = Nv(2)*(i-0)+j;   % > SE.
                            struct.ConnectivityList(Nc(2)*(i-1)+j,3) = Nv(2)*(i-0)+j+1; % > NE.
                            struct.ConnectivityList(Nc(2)*(i-1)+j,4) = Nv(2)*(i-1)+j+1; % > NW.
                        end
                    end
                    %  > Points.
                    struct.Points = cat(2,reshape(Xd,[],1),reshape(Yd,[],1));
                case "v"
                    % >> w/ triangles.
                    struct = delaunayTriangulation(reshape(Xd,[],1),reshape(Yd,[],1));
                otherwise
                    return;
            end
        end
        %  > 3.1.1. -------------------------------------------------------
        function [Xd,Yd] = T1(Xt,Yt)
            %  > f.
            f {1} = @(u,v) u;
            f {2} = @(u,v) v;
            %  > X/Yd.
            Xd    = f{1}(Xt,Yt);
            Yd    = f{2}(Xt,Yt);
        end
        %  > 3.1.2. -------------------------------------------------------
        function [Xd,Yd] = T2(Xt,Yt)
            %  > f.
            c (1) = 2.*pi;
            c (2) = 1./(2.*c(1));
            f {1} = @(u,v) u+c(2).*sin(c(1).*u).*sin(c(1).*v);
            f {2} = @(u,v) v+c(2).*sin(c(1).*u).*sin(c(1).*v);
            %  > X/Yd.
            Xd    = f{1}(Xt,Yt);
            Yd    = f{2}(Xt,Yt);
        end
        %  > 3.1.3. -------------------------------------------------------
        function [Xd,Yd] = T3(Xt,Yt)
            %  > f.
            c (1) = 0.75;
            f {1} = @(u,v) u;
            f {2} = @(u,v) v-c(1)/pi*sin(pi*v);
            %  > X/Yd.
            Xd    = f{1}(Xt,Yt);
            Yd    = f{2}(Xt,Yt);
        end
        %  > 3.1.4. -------------------------------------------------------
        function [Xd,Yd] = T4(Xt,Yt)
            %  > Auxiliary variables.
            [X(1),X(2)] = MinMaxElem(Xt); L = X(2)-X(1); Nv(1) = size(Xt,2); Nc(1) = Nv(1)-1;
            [Y(1),Y(2)] = MinMaxElem(Yt); H = Y(2)-Y(1); Nv(2) = size(Yt,1); Nc(2) = Nv(2)-1;
            %  > f.
            c (1)   = 0.1;
            c (2)   = 1.0.*pi;
            f {1}   = @(i,j) L.*(i-1)./Nc(1).*(c(1).*(Nv(1)-i)./Nc(1).*sin(c(2).*H./L.*(j-1)./Nc(2))+1);
            f {2}   = @(i,j) H.*(j-1)./Nc(2).*(c(1).*(Nv(2)-j)./Nc(2).*sin(c(2).*H./L.*(i-1)./Nc(1))+1);
            %  > X/Yd.
            [it,jt] = meshgrid(1:Nv(1),1:Nv(2));
            Xd      = f{1}(it,jt);
            Yd      = f{2}(it,jt);
        end
        
        %% > 4. -----------------------------------------------------------
        %  > Analytic functions.
        % >> 4.1. ---------------------------------------------------------
        function [fh] = func(t,v)
            %  > Symbolic variables.
            syms x y;
            
            switch t
                case 1
                    xc   = v(1);
                    yc   = v(2);
                    i    = v(3);
                    func = exp(-i.*((x-xc).^2+(y-yc).^2));
                otherwise
                    return;
            end
            fh = matlabFunction(func,'Vars',{x,y});
        end
    end
end