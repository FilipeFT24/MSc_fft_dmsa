classdef B_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Compute analytic functions/values/source term.
        function [obj_f] = Update_func(inp,msh)
            %  > Auxiliary variables.
            Xc = msh.c.Xc;
            Nc = msh.c.Nc;
            Xv = msh.f.Xv;
            v  = inp.pv.v(1);
            g  = inp.pv.v(2);
            
            % >> Set...
            %  > ...analytic function.
            syms x;
            f_1      = inp.pv.f(x);
            f_2      = diff(f_1,x);
            f_3      = diff(f_2,x);
            h        = v.*f_2-g.*f_3;
            %  > ... function handles.
            f  {1}   = matlabFunction(f_1);
            f  {2}   = matlabFunction(f_2);
            i        = matlabFunction(int(h));
            % >> Compute...
            %  > ...analytic cell/face values.
            a.c(:,1) = f{1}(Xc);
            a.f(:,1) = f{1}(Xv);
            a.f(:,2) = f{2}(Xv);
            %  > ...analytic (volumetric) source term.
            j        = 1:Nc;
            fv (j,1) = i(Xv(j+1))-i(Xv(j));

            % >> Update...
            obj_f.av = a;
            obj_f.fh = f;
            obj_f.fv = fv;
            obj_f.ps = B_1_1D.pol_shape;
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Create function handle of \phi and \nabla\phi with 100 terms (polynomial shape). 
        function [func] = pol_shape()
            syms x f;
            n       = 100;
            l_1     = 1:n;
            func{1} = (x-f).^(l_1-1);
            func{1} = matlabFunction(func{1});
            l_2     = 1:n-1;
            func{2} = l_2.*(x-f).^(l_2-1);
            func{2} = matlabFunction(func{2});
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