classdef B_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        %  > Compute analytic functions/values/source term.
        function [pde] = Update_pde(inp,msh)
            %  > Auxiliary variables.
            v  = inp.pv.v(1);
            g  = inp.pv.v(2);
            
            %  > Set analytic function.
            syms x;
            f   {1} = inp.pv.f(x);
            f   {2} = diff(f{1},x);
            func    = v.*f{2}-g.*diff(f{2},x);
            fn.f{1} = matlabFunction(f{1});
            fn.f{2} = matlabFunction(f{2});
            fn.i    = matlabFunction(int(func));
            
            % >> Compute analytic...
            %  > ...cell values.
            av.c  (:,1) = fn.f{1}(msh.c.Xc);
            %  > ...face values.
            av.f  (:,1) = fn.f{1}(msh.f.Xv);
            av.f  (:,2) = fn.f{2}(msh.f.Xv);
            %  > ...volumetric source term.
            i           = 1:msh.c.Nc;
            A     (i)   = msh.f.Xv(i);
            B     (i)   = msh.f.Xv(i+1);
            fn.vol(i,1) = fn.i(B(i))-fn.i(A(i));
            
            %  > Update 'pde' structure.
            pde.av = av;
            pde.fn = fn;
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