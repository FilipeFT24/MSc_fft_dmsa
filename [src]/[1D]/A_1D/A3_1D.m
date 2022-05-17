classdef A3_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Update field 'f.fh' (function handles).
        function [fh] = Update_fh(inp)
            %  > Auxiliary variables.
            syms x;
            
            % >> Analytic function/derivatives.
            %  > Function.
            f_d{1} = diff(inp.f(x),x);
            f_d{2} = diff(f_d{1},x);
            func   = inp.c(1).*f_d{1}+inp.c(2).*f_d{2};
            intc   = int(func);
            %  > Function derivatives up to (+1(boundary)+2/4(next stencil)=+3/5)...
            if inp.t_terms.allow && ~inp.p_adapt.allow
                m = length(inp.c);
                n = inp.t_terms.n;
                i = 1:m;
                j = inp.p(i)+n(i);
                for k = 1:max(j)+5
                    d_f{k} = diff(inp.f(x),x,k)/factorial(k);
                    d_f{k} = matlabFunction(d_f{k});
                end
                fh.f.d2    = d_f;
            end
            %  > f.
            fh.f.d1{1} = matlabFunction(f_d{1});
            fh.f.d1{2} = matlabFunction(f_d{2});
            fh.f.f     = matlabFunction(inp.f(x));
            %  > func.
            fh.func.f = matlabFunction(func);
            fh.func.i = matlabFunction(intc);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Update field 'f.av' (cell/face values).
        function [av] = Update_av(msh,f)
            %  > Auxiliary variables.
            Nc = msh.c.Nc;
            Xc = msh.c.Xc;
            Xv = msh.f.Xv;

            av.c(:,1) = f.f    (Xc);
            av.f(:,1) = f.f    (Xv);
            av.f(:,2) = f.d1{1}(Xv);
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Update field 'f.bd' (boundary face values).
        function [bd] = Update_bd(inp,msh,f)
            %  > Auxiliary variables.
            Nf    = msh.f.Nf;
            Xv    = msh.f.Xv;
            
            bd.i  = [1,Nf];
            bd.t  = inp.b.t;
            bd.x  = [Xv(1),Xv(end)];
            for i = 1:length(bd.t)
                switch bd.t(i)
                    case "Dirichlet"
                        bd.v(i) = f.f    (bd.x(i));
                    case "Neumann"
                        bd.v(i) = f.d1{1}(bd.x(i));
                    case "Robin"
                        bd.v(i) = f.f    (bd.x(i))-inp.c(2)./inp.c(1).*f.d1{1}(bd.x(i));
                    otherwise
                        return;
                end
            end
        end
        % >> 1.4. ---------------------------------------------------------
        %  > Update field 'f.st' ((volumetric) source term).
        function [st] = Update_st(msh,f)
            %  > Auxiliary variables.
            Nc = msh.c.Nc;
            Xv = msh.f.Xv;
            
            i             = 1:Nc;
            st(i,1)       = f(Xv(i+1))-f(Xv(i));
            st(isinf(st)) = 0;
            st(isnan(st)) = 0;
        end        
    end
end