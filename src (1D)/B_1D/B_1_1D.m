classdef B_1_1D
    methods (Static)
        %% > Wrap-up B_1 (1D).
        function [pde] = WrapUp_B_1_1D(msh,v,g,np,p_adapt,ft)
            % >> Compute...
            %  > ...analytic functions/values.
            pde.fn = B_1_1D.Set_fn(msh,v,g,ft);
            pde.sn = B_1_1D.Compute_blkf_blkc(msh,pde.fn);
            %  > ...source term.
            pde.FV = B_1_1D.Compute_ST(msh,np,p_adapt,pde.fn);
        end
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [fn] = Set_fn(msh,v,g,ft)
            switch char(ft)
                case 'sin'
                    syms x;
                    i    = 9;
                    f{1} = sin(i.*pi.*x);
                case 'exp'
                    syms x;
                    c    = 1./2.*(max(msh.f.Xv)-min(msh.f.Xv));
                    i    = 100;
                    f{1} = exp(-i.*((x-c).^2));
                otherwise
                    return;
            end
            f{2} = diff(f{1},x);
            func = v.*f{2}-g.*diff(f{2},x);
            intf = int(func);
            
            % >> Set...
            %  > ...symbolic functions/handles.
            fn.sym{1} = f{1};
            fn.sym{2} = f{2};
            fn.f  {1} = matlabFunction(f{1});
            fn.f  {2} = matlabFunction(f{2});
            fn.func   = matlabFunction(func);
            fn.int    = matlabFunction(intf);
        end
        % >> 1.2. ---------------------------------------------------------
        function [s] = Compute_blkf_blkc(msh,fn)
            %  > Cell values.
            i        = 1:msh.c.NC;
            s.c(i,1) = fn.f{1}(msh.c.Xc(i));
            s.c(i,2) = fn.f{2}(msh.c.Xc(i));
            
            %  > Face values.
            j        = 1:msh.f.NF;
            s.f(j,1) = fn.f{1}(msh.f.Xv(j));
            s.f(j,2) = fn.f{2}(msh.f.Xv(j));
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [x,j,Q_1D] = CD_1D(ng)
            %  > x(csi) = a*(1-csi)/2+b*(1+csi)/2.
            %  > j(csi) = d(x)/d(csi) = (b-a)/2.
            syms a b csi;
            x    = a.*(1-csi)./2+b.*(1+csi)./2;
            x    = matlabFunction(x);
            j    = (b-a)./2;
            j    = matlabFunction(j);
            %  > Quadrature abcissas/weights.
            Q_1D = quadGaussLegendre(ng);
        end
        % >> 2.2. ---------------------------------------------------------
        function [I] = ApproxIntegral(a,b,f,x,j,Q_1D)
            I = 0;
            for i = 1:length(Q_1D.Points)
                x_CD(i) = x(a,b,Q_1D.Points(i));
                I       = I+Q_1D.Weights(i).*j(a,b).*f(x_CD(i));
            end
        end
        % >> 2.3 ----------------------------------------------------------
        function [FV] = Compute_ST(msh,np,p_adapt,fn)
            %  > Interval extrema.
            i    = 1:msh.c.NC;
            A(i) = msh.f.Xv(i);
            B(i) = msh.f.Xv(i+1);
            
            if ~p_adapt
                %  > Appoximated integral.
                [x,j,Q_1D] = B_1_1D.CD_1D(np);
                func  = fn.func;
                for i = 1:msh.c.NC
                    FV(i,1) = B_1_1D.ApproxIntegral(A(i),B(i),func,x,j,Q_1D);
                end
            else
                %  > Exact integral.
                func  = fn.int;
                for i = 1:msh.c.NC
                    FV(i,1) = func(B(i))-func(A(i));
                end
            end
        end
    end
end