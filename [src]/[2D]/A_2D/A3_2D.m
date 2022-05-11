classdef A3_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "e" (error).
        function [f] = Initialize_f(inp,msh)
            %  > "fh" (function handles).
            f.fh = A3_2D.Update_fh(inp);
            %  > "av" (analytic values).
            f.av = A3_2D.Update_av(msh,f.fh.f);
            %  > "bd" (boundary values).
            f.bd = A3_2D.Update_bd(inp,msh,f.fh.f);
            %  > "st" (source term).
            f.st = A3_2D.Update_st(msh,f.fh.func.f);
            %  > "qd" (1D quadrature).
            for i = 1:size(inp.p.p,1)
                for j = 1:size(inp.p.p,2)
                    f.qd{i,j} = A3_2D.Q_1D(inp.p.p(i,j));
                end
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update field "f.fh" (function handles).
        function [fh] = Update_fh(inp)
            %  > Symbolic variables.
            x = sym  ('x');
            y = sym  ('y');
            v = inp.c(1,:);
            g = inp.c(2,:);
            
            %  > f.
            f            = inp.f([x,y]);
            df       (1) = diff (f,x);
            fh.f.d   {1} = matlabFunction(df(1),'Vars',{[x,y]});
            df       (2) = diff (f,y);
            fh.f.d   {2} = matlabFunction(df(2),'Vars',{[x,y]});
            fh.f.f       = matlabFunction(f    ,'Vars',{[x,y]});
            %  > func.
            func         = v{1}([x,y]).*df(1)+v{2}([x,y]).*df(2)+g{1}([x,y]).*diff(df(1),x)+g{2}([x,y]).*diff(df(2),y);
            fh.func.f    = func;
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update field "f.av" (analytic cell/face values).
        function [av] = Update_av(msh,f)
            av.c   (:,1) = f.f   (msh.c.c.xy.c);
            av.f{1}(:,1) = f.f   (msh.f.xy.c);
            av.f{2}(:,1) = f.d{1}(msh.f.xy.c);
            av.f{2}(:,2) = f.d{2}(msh.f.xy.c);
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Update field "f.bd" (boundary face indices/values).
        function [bd] = Update_bd(inp,msh,f)
            %  > (Boundary) face indices.
            bd.i(:,1) = find(~msh.f.logical);
            %  > Identify...
            for i = 1:numel(bd.i)
                c       = msh.f.ic  {bd.i(i,1)};
                Sf(i,:) = msh.c.f.Sf{c}(msh.c.f.if(c,:) == bd.i(i,1),:);
                %  > Select...
                if     Sf(i,1) >  0 && Sf(i,2) == 0, t(i,1) = 1; %  > East.
                elseif Sf(i,1) == 0 && Sf(i,2) >  0, t(i,1) = 2; %  > North.
                elseif Sf(i,1) <  0 && Sf(i,2) == 0, t(i,1) = 3; %  > West.
                elseif Sf(i,1) == 0 && Sf(i,2) <  0, t(i,1) = 4; %  > South.
                else
                    return;
                end
            end
            for i = 1:numel(bd.i)
                bd_j{i} = find(t == i);
            end
            %  > (Boundary) face values.
            for i = 1:numel(inp.b.t)
                bd_fi = bd_j{i};
                xy_fc = msh.f.xy.c(bd.i(bd_j{i}),:);
                switch inp.b.t(i)
                    case "Dirichlet"
                        %  > V = \phi(f).
                        bd.v(bd_fi,1) = f.f(xy_fc);
                    case "Neumann"
                        %  > V = (\nabla\phi(f)_x,\nabla\phi(f)_y)*Sf.
                        for j  = 1:numel(bd_fi)
                            bd.v(bd_fi(j),1) = ...
                                [f.d{1}(xy_fc(j,:)),f.d{2}(xy_fc(j,:))]*Sf(bd_fi(j),:)';
                        end
                    case "Robin"
                        %  > Evaluate convective/diffusive coefficients at face centroid...
                        for j = 1:size(inp.c,1)
                            for k = 1:size(inp.c,2)
                                c_fc{j}(:,k) = inp.c{j,k}(xy_fc);
                            end
                        end
                        GoV = c_fc{2}./c_fc{1};
                        %  > V = [GoV.*(\nabla\phi(f)_x,\nabla\phi(f)_y)]*Sf-\phi(f)*Sf.
                        for j  = 1:numel(bd_fi)
                            bd.v(bd_fi(j),1) = ...
                                [f.f(xy_fc(j,:)),f.f(xy_fc(j,:))]*Sf(bd_fi(j),:)'-GoV(j,:).*[f.d{1}(xy_fc(j,:)),f.d{2}(xy_fc(j,:))]*Sf(bd_fi(j),:)';
                        end
                end
            end
        end
        % >> 1.5. ---------------------------------------------------------
        %  > Update field "f.st" ((volumetric) source term).
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
            st(isinf(st)) = 0;
            st(isnan(st)) = 0;
        end
        %  > 1.5.1. -------------------------------------------------------
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
        %  > 1.5.2. -------------------------------------------------------
        %  > Auxiliary function #2: set up integrand function (int_c).
        function [int_c] = int(c,z,f)
            %  > Substitute...
            sub_c = subs(f,{sym('x'),sym('y')},{z(1),z(2)});
            %  > Convert to function handle...
            int_c = matlabFunction(sub_c,'Vars',{c,sym('u'),sym('v')});
            int_c = @(c)(@(u,v) int_c(c,u,v));
        end
        %  > 1.5.3. -------------------------------------------------------
        %  > Auxiliary function #3: compute (cell) jacobian (jac_c) and its absolute determinat (d_abs).
        function [d_abs] = jac(c,z)
            %  > Jacobian.
            jac_c = jacobian(z,[sym('u'),sym('v')]);
            %  > Jacobian (absolute) determinant.
            d_abs = abs(det(jac_c));
            d_abs = matlabFunction(d_abs,'Vars',{c});
        end
        %  > 1.5.4. -------------------------------------------------------
        %  > Auxiliary function #4: compute (cell) double integral ("v").
        function [st] = st_aux(map,xy_cv)
            st = map.d(xy_cv).*integral2(map.i(xy_cv),0,1,0,@(u) 1-u);
        end
        %  > 1.5.5. -------------------------------------------------------
        %  > Auxiliary function #5: 1D quadrature: x(u) = x(1)*(1-u)/2+x(2)*(1+u)/2.
        %                                          j(u) = d(x)/d(csi).
        function [qd] = Q_1D(p)
            %  > Symbolic variables (u,x).
            u = sym('u');
            x = sym('x',[1,2]);
            
            %  > x(u) and j(u).
            qd.xu = @(u,x) x(1).*(1-u)./2+x(2).*(1+u)./2;
            ju    = jacobian(qd.xu(u,x),u);
            qd.ju = matlabFunction(ju,'Vars',{x});
            %  > Quadrature points/weights (Gauss-Legendre).
            qd.pw = quadGaussLegendre(ceil(p./2));
        end
    end
end