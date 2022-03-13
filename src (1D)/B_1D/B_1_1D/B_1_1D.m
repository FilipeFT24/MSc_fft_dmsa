classdef B_1_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        %  > Compute analytic functions/values/source term.
        function [pde] = Update_pde(inp,msh)
            %  > Auxiliary variables.
            v = inp.pv.v(1);
            g = inp.pv.v(2);
            
            %  > Set analytic function.
            syms x;
            switch char(inp.pv.f)
                case "1"
                    i    = 6;
                    f{1} = sin(i.*pi.*x);
                case "2"
                    c    = 1./2.*(max(msh.f.Xv)-min(msh.f.Xv));
                    i    = 50;
                    f{1} = exp(-i.*((x-c).^2));
                otherwise
                    return;
            end
            f   {2} = diff(f{1},x);
            func    = v.*f{2}-g.*diff(f{2},x);
            fn.f{1} = matlabFunction(f{1});
            fn.f{2} = matlabFunction(f{2});
            fn.i    = matlabFunction(int(func));
            
            % >> Compute analytic...
            %  > ...cell values.
            av.c (:,1) = fn.f{1}(msh.c.Xc);
            %  > ...face values.
            av.f (:,1) = fn.f{1}(msh.f.Xv);
            av.f (:,2) = fn.f{2}(msh.f.Xv);
            %  > ...source term.
            i          = 1:msh.c.NC;
            A    (i)   = msh.f.Xv(i);
            B    (i)   = msh.f.Xv(i+1);
            fn.st(i,1) = fn.i(B(i))-fn.i(A(i));
            
            %  > Update/sort 'pde' structure.
            pde.av = av;
            pde.fn = fn;
        end
        
        %% > 2. -----------------------------------------------------------
        %  > 1D GQ functions (NOT USED).
        % >> 2.1. ---------------------------------------------------------
        function [x,j,Q_1D] = GQ_1D(ng)
            syms a b csi;
            x    = a.*(1-csi)./2+b.*(1+csi)./2;
            x    = matlabFunction(x);
            j    = (b-a)./2;
            j    = matlabFunction(j);
            Q_1D = quadGaussLegendre(ng);
        end
        % >> 2.2. ---------------------------------------------------------
        function [I] = Ap_I(a,b,f,x,j,Q_1D)
            I = 0;
            for i = 1:length(Q_1D.Points)
                x_CD(i) = x(a,b,Q_1D.Points(i));
                I       = I+Q_1D.Weights(i).*j(a,b).*f(x_CD(i));
            end
        end
    end
end