classdef A3_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Initialize field "f" (analytic function(s), boundary values, etc.).
        function [f] = Initialize_f(inp,msh)
            %  > "fh" (function handles).
            f.fh = A3_1D.Update_fh(inp);
            %  > "bd" (boundary values).
            f.bd = A3_1D.Update_bd(inp,msh,f.fh.f);
            %  > "st" (source term).
            f.st = A3_1D.Update_st(msh,f.fh.func.i);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update field "fh" (function handles).
        function [fh] = Update_fh(inp)
            %  > Symbolic variables.
            x = sym('x');
            c = inp.c;
            
            %  > "f".
            f            = inp.f(x);
            df       {1} = diff (f,x);
            df       {2} = diff (df(1),x);
            fh.f.d       = matlabFunction(df{1});
            fh.f.f       = matlabFunction(f);
            %  > "func".
            func         = c{1}.*df{1}+c{2}.*df{2};
            intc         = int(func);
            fh.func.f    = matlabFunction(func);
            fh.func.i    = matlabFunction(intc);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update field "bd" (boundary face indices/values).
        function [bd] = Update_bd(inp,msh,f)
            %  > "i", "t" and "x".
            bd.i  = [1,msh.f.Nf];
            bd.t  = inp.b.t;
            bd.x  = msh.f.Xv(bd.i)';
            %  > "v".
            for i = 1:numel(bd.t)
                switch bd.t(i)
                    case "Dirichlet"
                        bd.v(i) = f.f(bd.x(i));
                    case "Neumann"
                        bd.v(i) = f.d(bd.x(i));
                    case "Robin"
                        bd.v(i) = f.f(bd.x(i))+inp.c{2}(bd.x(i))./inp.c{1}(bd.x(i)).*f.d(bd.x(i));
                    otherwise
                        return;
                end
            end
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Update field "st" ((volumetric) source term).
        function [st] = Update_st(msh,f)
            for i = 1:msh.c.Nc
                st(i,1) = f(msh.f.Xv(i+1))-f(msh.f.Xv(i));
            end
            st(isinf(st) | isnan(st)) = 0;
        end
    end
end