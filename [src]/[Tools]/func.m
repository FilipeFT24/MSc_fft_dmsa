classdef func
    methods (Static)
        %% > 1. -----------------------------------------------------------
        %  > Set "struct" structure (msh).
        % >> 1.1. ---------------------------------------------------------
        function [struct] = Set_struct(inp)
            %  > Limits.
            for i = 1:size(inp.Lim,1)
                [L(1,i),L(2,i)] = MinMaxElem(inp.Lim(i,:));
            end
            %  > Select cell polyhedral (type)...
            switch inp.p
                % >> w/ squares.
                case "s"
                    %  > XYv.
                    for i = 1:size(inp.Lim,1)
                        Nv (i) = round(diff(L(:,i))./inp.h)+1; 
                        Nc (i) = Nv(i)-1;
                        XYv{i} = linspace(L(1,i),L(2,i),Nv(i));
                    end
                    %  > XYt,d.
                    [Xt,Yt] = meshgrid(XYv{1},XYv{2});
                    switch inp.t
                        case 0, [Xd,Yd] = func.demo_0(Xt,Yt);
                        case 1, [Xd,Yd] = func.demo_1(Xt,Yt);
                        case 2, [Xd,Yd] = func.demo_2(Xt,Yt,inp.h);
                        otherwise
                            return;
                    end
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
                % >> w/ triangles.
                case "v"
                    [struct.Points,struct.ConnectivityList] = ...
                        distmesh2d(@(p) drectangle(p,L(1,1),L(2,1),L(1,2),L(2,2)),@(p) ones(size(p,1),1),inp.h,L,[L(1,:);diag(L)';diag(fliplr(L'),0)';L(2,:)],[1.0e-3,1.0e-3]);
                otherwise
                    return;
            end
        end
        %  > 1.1.1. -------------------------------------------------------
        function [Xd,Yd] = demo_0(Xt,Yt)
            %  > f.
            f {1} = @(u,v) u;
            f {2} = @(u,v) v;
            %  > X/Yd.
            Xd    = f{1}(Xt,Yt);
            Yd    = f{2}(Xt,Yt);
        end
        %  > 1.1.2. -------------------------------------------------------
        function [Xd,Yd] = demo_1(Xt,Yt)
            %  > Auxiliary variables.
            [X(1),X(2)] = MinMaxElem(Xt); L = X(2)-X(1); Nv(1) = size(Xt,2); Nc(1) = Nv(1)-1;
            [Y(1),Y(2)] = MinMaxElem(Yt); H = Y(2)-Y(1); Nv(2) = size(Yt,1); Nc(2) = Nv(2)-1;
            %  > f.
            c (1)   = 0.1;
            c (2)   = 10.*pi;
            f {1}   = @(i,j) L.*(i-1)./Nc(1).*(c(1).*(Nv(1)-i)./Nc(1).*sin(c(2).*H./L.*(j-1)./Nc(2))+1);
            f {2}   = @(i,j) H.*(j-1)./Nc(2).*(c(1).*(Nv(2)-j)./Nc(2).*sin(c(2).*H./L.*(i-1)./Nc(1))+1);
            %  > X/Yd.
            [it,jt] = meshgrid(1:Nv(1),1:Nv(2));
            Xd      = f{1}(it,jt);
            Yd      = f{2}(it,jt);
        end
        function [Xd,Yd] = demo_2(Xt,Yt,h)
            %  > Auxiliary variables.
            k       = 0.1;
            m       = size(Xt,1)-2; i = 2:m+1;
            n       = size(Xt,2)-2; j = 2:n+1;
            %  > X/Yd.
            Xt(i,j) = Xt(i,j)+randn(m,n).*k.*h; Xd = Xt;
            Yt(i,j) = Yt(i,j)+randn(m,n).*k.*h; Yd = Yt;
        end
               
        %% > 2. -----------------------------------------------------------
        % >> Sort structures "msh" (1D/2D).
        % >> 2.1. ---------------------------------------------------------
        function [msh] = Sort_1D_msh(msh)
            % >> msh.
            msh     = orderfields(msh    ,{'c','d','f'});
            %  > c.
            msh.c   = orderfields(msh.c  ,{'f','Nc','Volume','Xc'});
            msh.c.f = orderfields(msh.c.f,{'if','Sf'});
            %  > d.
            msh.d   = orderfields(msh.d  ,{'h'});
            %  > f.
            msh.f   = orderfields(msh.f  ,{'ic','Nf','Xv'});
        end
        % >> 2.2. ---------------------------------------------------------
        function [msh] = Sort_2D_msh(msh)
            % >> msh.
            msh        = orderfields(msh       ,{'c','d','f','v','struct'});
            %  > c.
            msh.c      = orderfields(msh.c     ,{'c','f','h','logical','Nc','Volume'});
            msh.c.c    = orderfields(msh.c.c   ,{'nb','xy'});
            msh.c.c.nb = orderfields(msh.c.c.nb,{'f','v'});
            msh.c.c.xy = orderfields(msh.c.c.xy,{'c','v'});
            msh.c.f    = orderfields(msh.c.f   ,{'if','Sf'});
            msh.c.h    = orderfields(msh.c.h   ,{'h','xy'});
            %  > d.
            msh.d      = orderfields(msh.d     ,{'h'});
            %  > f.
            msh.f      = orderfields(msh.f     ,{'ic','iv','logical','Nf','xy'});
            msh.f.xy   = orderfields(msh.f.xy  ,{'c','v'});
            %  > v.
            msh.v      = orderfields(msh.v     ,{'ic','if','logical'});
        end
        
        %% > 3. -----------------------------------------------------------
        % >> Other functions.
        % >> 3.1. ---------------------------------------------------------
        %  > Similar (faster) version of other built-in functions (to compute distance between 2 points).
        function [D] = dist(XY)
            D = sqrt((XY(1,1)-XY(2,1)).^2+(XY(1,2)-XY(2,2)).^2);
        end
        % >> 3.2. ---------------------------------------------------------
        %  > Similar (faster) version of built-in function "mean".
        function [y] = mean(x,j)
            y = sum(x,j,'default','includenan')./size(x,j);
        end
        % >> 3.3. ---------------------------------------------------------
        %  > Similar (faster) version of built-in function "setdiff" (w/o verification(s)).
        function [z] = setdiff(x,y)
            w    = false(1,max(max(x),max(y)));
            w(x) = true;
            w(y) = false;
            z    = x(w(x));
        end
        % >> 3.4. ---------------------------------------------------------
        %  > Similar (faster) version of built-in function "meshgrid" (w/o verification(s)).
        function [a,b] = meshgrid(v)
            b = v'*ones(1,numel(v));
            a = b';
        end
        % >> 3.5. ---------------------------------------------------------
        %  > Similar (faster) version of following expression: arrayfun(@(x) find(a == x),b);
        function [j] = find_c(a,b)
            i = repmat(a,1,numel(b)) == b'; [j,~] = find(i);
        end
        % >> 3.6. ---------------------------------------------------------
        %  > 3.6.1. -------------------------------------------------------
        %  > Rescaling to Solve a Linear System (Ax=b).
        %  > Reference: https://www.mathworks.com/help/matlab/ref/equilibrate.html#mw_beaf3b94-7114-4652-b483-a36c2acf7b94
        function [x] = backslash(A,b)
            %  > Use built-in function "equilibrate".
            [P,R,C] = equilibrate(A);
            %  > Solve linear system...
            B = R*P*A*C;
            d = R*P*b;
            x = C*(B\d);
        end
        %  > 3.6.2. -------------------------------------------------------
        %  > Compute CLS matrices.
        %  > x = (H'H)^{-1}*(H'y-C'*(C*(H'H)^{-1}*C')^{-1}*(C*(H'*H)^{-1}H'y-b)).
        function [t] = cls_t(b,C,HTH,HT)
            %  > Auxiliary variables.
            C_HTH_CT    = C *func.backslash(HTH,C');
            C_HTH_HT    = C *func.backslash(HTH,HT);
            X           = C'*func.backslash(C_HTH_CT,C_HTH_HT);
            Y           = C'*func.backslash(C_HTH_CT,b);
            %  > CLS terms (cell-dependent matrix/bd_v-dependent vector).
            t       {1} =    func.backslash(HTH,HT-X);
            t       {2} =    func.backslash(HTH,Y);
        end
        % >> 3.7. ---------------------------------------------------------
        %  > Compute error norms (cell/face L1,L2 and L_infinity norms).
        function [L] = Set_n(E,V)
            if nargin == 1
                L(1,:) = func.mean(E,1);
                L(2,:) = func.mean(sqrt(E.^2),1);
                L(3,:) = max(E);
            else
                L(1,:) = sum(E.*V)./sum(V);
                L(2,:) = sum(sqrt((E.*V).^2))./sum(sqrt(V.^2));
                L(3,:) = max(E);
            end
        end
        % >> 3.8. ---------------------------------------------------------
        %  > Compute error slope.
        function [s] = Slope(h,e)
            [m,n]  = size(e);
            i      = 1:n;
            j      = 1:m-1;
            s(j,i) = log(e(j+1,i)./e(j,i))./log(h(j+1)./h(j));
        end
        % >> 3.9. ---------------------------------------------------------
        %  > Save .mat file.
        function [] = Save_mat(wd,V)
            save(strjoin([wd,".mat"],''),"V");
        end
    end
end