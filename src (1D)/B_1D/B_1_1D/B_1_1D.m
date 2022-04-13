classdef B_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Set analytic function/boundary values/(volumetric) source term/etc.
        function [obj_f] = Update_func(inp,msh)
            %  > Auxiliary variables.
            Nc = msh.c.Nc;
            Nf = msh.f.Nf;
            Xc = msh.c.Xc;
            Xv = msh.f.Xv;
            v  = inp.pv.v(1);
            g  = inp.pv.v(2);
            
            % >> Set...
            %  > ...analytic function handles.
            syms x;
            f_1      = inp.pv.f(x);
            f_2      = diff(f_1,x);
            f_3      = diff(f_2,x);
            h        = v.*f_2-g.*f_3;
            f    {1} = matlabFunction(f_1);
            f    {2} = matlabFunction(f_2);
            i        = matlabFunction(int(h));
            %  > ...analytic cell/face values.
            a.c(:,1) = f{1}(Xc);
            a.f(:,1) = f{1}(Xv);
            a.f(:,2) = f{2}(Xv);
            %  > ...analytic (volumetric) source term.
            j        = 1:Nc;
            fv (j,1) = i(Xv(j+1))-i(Xv(j));
            %  > ...boundary types/values.
            bd.t     = inp.pv.b;
            bd.x     = [Xv(1),Xv(Nf)];
            bd.v     = B_1_1D.Set_bd_v(bd.t,bd.x,f,g./v);
            
            % >> Update...
            obj_f.av = a;
            obj_f.fh = f;
            obj_f.fv = fv;
            obj_f.bd = bd;
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Set boundary values.
        function [bd_v] = Set_bd_v(bd_t,bd_x,f,gv)
            for i = 1:length(bd_t) 
                switch bd_t(i)
                    case "Dirichlet"
                        bd_v(i) = f{1}(bd_x(i));
                    case "Neumann"
                        bd_v(i) = f{2}(bd_x(i));
                    case "Robin"
                        bd_v(i) = f{1}(bd_x(i))+gv.*f{2}(bd_x(i));
                end
            end
        end
               
        %% > 2. -----------------------------------------------------------
        %  > 1D GQ functions (NOT USED).
        % >> 2.1. ---------------------------------------------------------
        function [gq] = GQ(n,f)
            syms a b csi;
            gq.n   = n;
            x      = a.*(1-csi)./2+b.*(1+csi)./2;
            gq.x   = matlabFunction(x);
            j      = (b-a)./2;
            gq.j   = matlabFunction(j);
            gq.GQ  = quadGaussLegendre(n);
            gq.cn  = factorial(n).^4./((2.*n+1).*(factorial(2.*n).^3));
            syms x;
            fn     = f;
            gq.fn  = matlabFunction(fn);
            dn     = diff(f,x,2.*n);
            gq.dn  = matlabFunction(dn);
        end
        % >> 2.2. ---------------------------------------------------------
        function [I] = Ap_I(gq,x)
            n = length(gq.GQ.Points);
            I = 0;
            for j = 1:n
                y(j) = gq.x(x(1),x(2),gq.GQ.Points(j));
                I    = I+gq.j(x(1),x(2)).*gq.GQ.Weights(j).*gq.fn(y(j));
            end
        end
        % >> 2.3. ---------------------------------------------------------
        function [en] = en(gq,x)
            en = gq.cn.*(x(2)-x(1)).^(2.*gq.n+1).*max(gq.dn(x(1)),gq.dn(x(2)));
        end
    end
end