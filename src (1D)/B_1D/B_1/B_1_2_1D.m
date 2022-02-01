classdef B_1_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1 ----------------------------------------------------------
        function [x,j,Q_1D] = CD_1D(ng)
            %  > Symbolic variables.
            syms a b csi;
            
            %  > x(csi) = a*(1-csi)/2+b*(1+csi)/2.
            x = a.*(1-csi)./2+b.*(1+csi)./2;
            x = matlabFunction(x);
            %  > j(csi) = d(x)/d(csi) = (b-a)/2.
            j = (b-a)./2;
            j = matlabFunction(j);
            %  > Quadrature abcissas/weights.
            Q_1D = quadGaussLegendre(ng);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1 ----------------------------------------------------------
        function [I] = ApproxIntegral(a,b,f,x,j,Q_1D)           
            I = 0;
            for i = 1:length(Q_1D.Points)
                x_CD(i) = x(a,b,Q_1D.Points(i));
                I       = I+Q_1D.Weights(i).*j(a,b).*f(x_CD(i));
            end
        end
        % >> 2.2 ----------------------------------------------------------
        function [F_Vol] = Compute_SourceTerm(msh,pde,st,ng)
            %  > Interval extrema.
            i    = 1:msh.c.NC;
            A(i) = msh.f.Xv(i);
            B(i) = msh.f.Xv(i+1);
            
            if ~st
                % >> Appoximated integral.
                %  > 1D Quadrature.
                [x,j,Q_1D] = B_1_2_1D.CD_1D(ng);
                %  > Source term.
                func  = pde.an.fn.func;
                for i = 1:msh.c.NC
                    F_Vol(i,1) = B_1_2_1D.ApproxIntegral(A(i),B(i),func,x,j,Q_1D);
                end
                F_Vol = sparse(F_Vol);
            else
                % >> Exact integral.
                %  > Source term.
                func  = pde.an.fn.int;
                for i = 1:msh.c.NC
                    F_Vol(i,1) = func(B(i))-func(A(i));
                end
                F_Vol = sparse(F_Vol);
            end
        end
    end
end