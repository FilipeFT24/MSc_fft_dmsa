classdef A3_2D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Update field 'f.fh' (function handles).
        function [fh] = Update_fh(inp)
            %  > Symbolic variables.
            syms x y;

            %  > f.
            f               = inp.f(x,y);
            df          (1) = diff (f,x);
            fh.f.d      {1} = matlabFunction(df(1),'Vars',{[x,y]});
            df          (2) = diff (f,y);
            fh.f.d      {2} = matlabFunction(df(2),'Vars',{[x,y]});
            fh.f.f          = matlabFunction(f    ,'Vars',{[x,y]});
            %  > func.
            v               = inp.c(1,:);
            g               = inp.c(2,:);
            func            = (v(1).*df(1)+v(2).*df(2))-(g(1).*diff(df(1),x)+g(2).*diff(df(2),y));
            fh.func.f       = func;
            fh.func.fh_f    = matlabFunction(func ,'Vars',{[x,y]});
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update field 'f.av' (cell/face values).
        function [av] = Update_av(msh,f)        
            av.c(:,1) = f.f   (msh.c.c.xy.c);
            av.f(:,1) = f.f   (msh.f.xy.c);
            av.f(:,2) = f.d{1}(msh.f.xy.c);
            av.f(:,3) = f.d{2}(msh.f.xy.c);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update field 'f.bd' (boundary face indices/values).
        function [bd] = Update_bd(inp,msh,f)
            bd.i(:,1) = find(~msh.f.logical);
            bd.v(:,1) = f.f ( msh.f.xy.c(~msh.f.logical,:)); 
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Update field 'f.st' ((volumetric) source term).
        function [st] = Update_st(msh,func,flag)
            %  > Auxiliary variables.
            map  = A3_2D.Map(func,flag);
            CL_c = msh.struct.ConnectivityList;

            for i = 1:msh.c.Nc
                %  > Split cell polyhedral into (#faces/vertices-2) triangles...
                for j = 1:size(msh.c.c.xy.v{i},1)-2
                    %  > Select vertex indices/coordinates.
                    xy_v  {i,1}{j} = msh.struct.Points(CL_c(i,[1,j+1,j+2]),:);
                    %  > Compute partial contributions.
                    st_aux{i,1}(j) = A3_2D.st_aux(map,xy_v{i,1}{j},flag);
                end
                st(i,1) = sum(st_aux{i,1}); 
            end
            st(isinf(st)) = 0;
            st(isnan(st)) = 0;
        end
        %  > 1.4.1. -------------------------------------------------------
        %  > Auxiliary function #1: 1) coordinate transformation ("v"): x(u,v) = a[1](1-u-v)+b[1](u)+c[1](v), where xv([1],[2],[3]) = (a[1],b[1],c[1]).
        %                                                               y(u,v) = a[2](1-u-v)+b[2](u)+c[2](v), where yv([1],[2],[3]) = (a[2],b[2],c[2]).
        function [map] = Map(func,flag)
            %  > Symbolic variables (u,v,c).
            syms u v;
            c = sym('c',[3,2]);
            
            %  > z = (x,y){u,v}.
            for k = 1:size(c,2)
                z(k) = c(1,k).*(1-u-v)+c(2,k).*u+c(3,k).*v;
            end
            %  > Assign to structure "map".
            if ~flag
                map.i = A3_2D.int(c,z,func);
            else
                ng    = 10;
                map.q = quadtriangle  (ng  ,'Domain',[0,0;1,0;0,1]);
                map.z = matlabFunction(z   ,'Vars'  ,{c,u,v});
                map.f = matlabFunction(func,'Vars'  ,{sym('x'),sym('y')});
            end
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
        function [st] = st_aux(map,xy_v,flag)
            if ~flag
                % >> w/ analytic integrand.
                st   = map.d(xy_v).*integral2(map.i(xy_v),0,1,0,@(u) 1-u);
            else
                % >> w/ 2D quadrature.
                %  > Points/weights/mapping...
                p    = map.q.Points;
                w    = map.q.Weights;
                z    = map.z(xy_v,p(:,1),p(:,2));
                %  > Individual contribution of each point/weight.
                st_i = map.f(z(:,1),z(:,2)).*w;
                %  > Global contribution.                
                st   = map.d(xy_v).*sum(st_i);
            end
        end
        %  > 1.4.4. -------------------------------------------------------
        %  > Auxiliary function #5: 1D Quadrature (1).
        function [Q_1D] = Q_1D(p)
            %  > x(u) = a*(1-u)/2+b*(1+u)/2.
            Q_1D.xu = @(x,u) x(1).*(1-u)./2+x(2).*(1+u)./2; 
            %  > j(u) = d(x)/d(csi) = (b-a)/2.
            Q_1D.j  = @(x)   (x(1)-x(2))./2;
            %  > Quadrature points/weights (Gauss-Legendre).
            Q_1D.pw = quadGaussLegendre(p); 
            Q_1D.pw.Weights = 1./2.*Q_1D.pw.Weights;
        end
    end
end