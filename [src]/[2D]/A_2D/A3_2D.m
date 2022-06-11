classdef A3_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "e" (error).
        function [f] = Initialize_f(inp,msh)
            %  > "fh" (function handles).
            f.fh = A3_2D.Update_fh(inp);
            %  > "bd" (boundary values).
            f.bd = A3_2D.Update_bd(inp,msh,f.fh.f);
            %  > "st" (source term).
            f.st = A3_2D.Update_st(msh,f.fh.func.f);
            %  > "qd" (1D quadrature).
            f.qd = A3_2D.Q_1D_1;
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update field "fh" (function handles).
        function [fh] = Update_fh(inp)
            %  > Symbolic variables.
            x = sym  ('x');
            y = sym  ('y');
            v = inp.c(1,:);
            g = inp.c(2,:);
            
            %  > "f".
            f              = inp.f([x,y]);
            df       {1,1} = diff (f,x);
            df       {1,2} = diff (f,y);
            df       {2,1} = diff (df{1,1},x);
            df       {2,2} = diff (df{1,2},y);
            fh.f.d     {1} = matlabFunction(df{1,1},'Vars',{[x,y]});
            fh.f.d     {2} = matlabFunction(df{1,2},'Vars',{[x,y]});
            fh.f.f         = matlabFunction(f      ,'Vars',{[x,y]});
            %  > "func".
            fh.func.f      = v{1}([x,y]).*df{1,1}+v{2}([x,y]).*df{1,2}+...
                g{1}([x,y]).*df{2,1}+g{2}([x,y]).*df{2,2};
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update field "bd" (boundary face indices/type).
        %  > NOTE: hard coded for square domain.
        function [bd] = Update_bd(inp,msh,fh)
            %  > "i".
            bd.i(:,1) = find(~msh.f.logical);
            %  > "t" and "v".
            for i  = 1:numel(bd.i)
                c  = msh.f.ic  {bd.i(i,1)};
                Sf = msh.c.f.Sf{c}(msh.c.f.if(c,:) == bd.i(i,1),:);
                %  > Identify boundary type...
                if     Sf(1) >  0 && Sf(2) == 0, bd.t(i,1) = 1; %  > East (E).
                elseif Sf(1) == 0 && Sf(2) >  0, bd.t(i,1) = 2; %  > North(N).
                elseif Sf(1) <  0 && Sf(2) == 0, bd.t(i,1) = 3; %  > West (W).
                elseif Sf(1) == 0 && Sf(2) <  0, bd.t(i,1) = 4; %  > South(S).
                else
                    return;
                end
                %  > Compute boundary value.
                xf_c = msh.f.xy.c(bd.i(i),:);
                switch inp.b.t(bd.t(i))
                    case "Dirichlet"
                        bd.v(i,1) = fh.f(xf_c);
                    case "Neumann"
                        bd.v(i,1) = Sf*[fh.d{1}(xf_c);fh.d{2}(xf_c)];
                    case "Robin"
                        for j = 1:size(inp.c,2)
                            GoV(j,1) = inp.c{2,j}(xf_c)./inp.c{1,j}(xf_c);
                        end
                        bd.v(i,1) = Sf*[fh.f(xf_c);fh.f(xf_c)]+Sf*[GoV(1).*fh.d{1}(xf_c);GoV(2).*fh.d{2}(xf_c)];
                end
            end
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Update field "st" ((volumetric) source term).
        function [st] = Update_st(msh,func)
            %  > Map...
            map = A3_2D.Map(func);
            
            for i = 1:msh.c.Nc
                for j = 1:size(msh.c.c.xy.v{i},1)-2
                    %  > Split cell polyhedral into #V-2 triangles and select vertex indices/coordinates...
                    xy_v  {i,1}{j} = msh.struct.Points(msh.struct.ConnectivityList(i,[1,j+1,j+2]),:);
                    %  > Compute partial contributions...
                    st_aux{i,1}(j) = A3_2D.st_aux(map,xy_v{i,1}{j});
                end
                st(i,1) = sum(st_aux{i,1});
            end
            st(isinf(st) | isnan(st)) = 0;
        end
        %  > 1.4.1. -------------------------------------------------------
        %  > Auxiliary function #1: (x,y)->(u,v): 1. x(u,v) = a[1](1-u-v)+b[1](u)+c[1](v), where xv([1],[2],[3]) = (a[1],b[1],c[1]).
        %                                         2. y(u,v) = a[2](1-u-v)+b[2](u)+c[2](v), where yv([1],[2],[3]) = (a[2],b[2],c[2]).
        function [map] = Map(f)
            %  > Symbolic variables (c,u,v).
            c = sym('c',[3,2]);
            u = sym('u');
            v = sym('v');
            
            %  > z = (x,y){u,v}.
            for k = 1:size(c,2)
                z(k) = c(1,k).*(1-u-v)+c(2,k).*u+c(3,k).*v;
            end
            %  > Assign to structure "map".
            map.i = A3_2D.int(c,z,f);
            map.d = A3_2D.jac(c,z);
        end
        %  > 1.4.2. -------------------------------------------------------
        %  > Auxiliary function #2: set up integrand function (int_c).
        function [int_c] = int(c,z,f)
            %  > Substitute...
            sub_c = subs(f,{sym('x'),sym('y')},{z(1),z(2)});
            %  > Convert to function handle...
            int_c = matlabFunction(sub_c,'Vars',{c,sym('u'),sym('v')});
            int_c = @(c)(@(u,v) int_c(c,u,v));
        end
        %  > 1.4.3. -------------------------------------------------------
        %  > Auxiliary function #3: compute (cell) jacobian (jac_c) and its absolute determinat (d_abs).
        function [d_abs] = jac(c,z)
            %  > Jacobian.
            jac_c = jacobian(z,[sym('u'),sym('v')]);
            %  > Jacobian (absolute) determinant.
            d_abs = abs(det(jac_c));
            d_abs = matlabFunction(d_abs,'Vars',{c});
        end
        %  > 1.4.4. -------------------------------------------------------
        %  > Auxiliary function #4: compute (cell) double integral ("v").
        function [st] = st_aux(map,xy_cv)
            st = map.d(xy_cv).*integral2(map.i(xy_cv),0,1,0,@(u) 1-u);
        end
        %  > 1.4.5. -------------------------------------------------------
        %  > 1.4.5.1. -----------------------------------------------------
        %  > Auxiliary function #5.1.
        function [qd] = Q_1D_1()
            qd.xu = @(u,x) x(1,:).*(1-u)./2+x(2,:).*(1+u)./2;
        end
        %  > 1.4.5.2. -----------------------------------------------------
        %  > Auxiliary function #5.2.
        function [Q] = Q_1D_2(n)
            %  > From "quadGaussLegendre(n,varargin)"...
            A       = zeros(1,n);
            B       = sqrt (((1:n-1)./(2:n))./((2*(0:n-2)+1)./((0:n-2)+1).*(2*(1:n-1)+1)./((1:n-1)+1)));
            J       = diag(B,1)+diag(A)+diag(B,-1);
            %  > Weights/points.
            [V,D]   = eig(J,'vector');
            [Q.x,I] = sort(D); Q.w = (2*V(1,I).^2)';
        end
    end
end