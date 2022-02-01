classdef SubClass_2_1
    methods (Static)
        %% > Wrap up SubClass_2_1.
        function [fn,blk_f,blk_c] = WrapUp_2_1(inp,msh)
            v = inp.pr.v;
            g = inp.pr.g;
            
            %  > #1: f,df,d2f analytic expression.
            fn = SubClass_2_1.Set_fn(v,g);
            %  > #2: Compute analytic solution.
            [blk_f,blk_c] = SubClass_2_1.Compute_blkf_blkc(msh,fn.f);
        end
        
        %% > Auxiliary functions.
        % >> Set face f, df and d2f.
        function [fn] = Set_fn(v,g)
            % >> Symbolic variable.
            syms x;
            
            %  > Phi.
            sig  = 0.10;
            mu   = 0.50;
            fn.f = 1./(sig.*sqrt(2.*pi)).*exp(-1./2.*((x-mu)./sig).^2);
            %  > grad(n)Phi.
            df   = diff(fn.f,x);
            d2f  = diff(diff(fn.f,x),x);
            %  > func.
            fn.func = v.*df-g.*d2f;
            
            %  > Function handles.
            [fn.f,fn.func] = deal(matlabFunction(fn.f),matlabFunction(fn.func));
        end
        
        % >> Evaluate at...
        function [blk_f,blk_c] = Compute_blkf_blkc(msh,f)
            %  > Face values.
            i          = 1:msh.f.NF;
            blk_f(i,1) = f(msh.f.Xv(i));
            %  > Cell values.
            j          = 1:msh.c.NC;
            blk_c(j,1) = f(msh.c.Xc(j));
        end
        
        %% > Tools.
        % >> #0.
        function [Pe] = Compute_Peclet
            Pe = abs(obj.V).*(obj.Xv_f-obj.Xv_i)./abs(obj.Gamma);
        end
        
        % >> #1.
        function [I] = ExactIntegral(f,a,b)
            %  > Symbolic variable(x).
            syms x;
            
            IX = int(f);
            I  = double(subs(IX,x,b))-double(subs(IX,x,a));
        end
        
        % >> #2.
        function [I] = ApproxIntegral(n,f,a,b)
            %  > Symbolic variable(x).
            syms x xi;
            
            x_xi = a.*(1-xi)./2+b.*(1+xi)./2;
            J    = (b-a)./2;
            
            %  > Location values/weights.
            [xg,wg] = lgwt(n,-1,1);
            for i = 1:size(xg,1)
                F_xi(i) = double(subs(x_xi,xi,xg(i)));
                I_i (i) = wg(i).*J.*f(F_xi(i));
            end
            I = sum(I_i);
        end
        
        % >> #3.
        function [F_Vol] = Compute_SourceTerm(n,fn,msh)
            for i = 1:msh.c.NC
                %  > Exact value (NOT USED).
                %F_Vol.Ex(i) = SubClass_2_1.ExactIntegral(fn.func,msh.f.Xv(i),msh.f.Xv(i+1));
                %  > Approximated value.
                F_Vol.Ap(i) = SubClass_2_1.ApproxIntegral(n,fn.func,msh.f.Xv(i),msh.f.Xv(i+1));
            end
            %F_Vol.Ex = sparse(F_Vol.Ex);
            F_Vol.Ap = sparse(F_Vol.Ap);
        end
    end
end