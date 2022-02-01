classdef SubClass_2_1
    methods (Static)
        %% > Wrap up SubClass_2_1.
        function [fn,bnd,blk] = WrapUp_2_1(msh)           
            %  > #1: f,df,d2f analytic expression.
            fn = SubClass_2_1.Set_fn(msh);
            %  > #2: Compute analytic solution.
            [bnd,blk] = SubClass_2_1.Compute_f_df_d2f(fn,msh);
        end
        
        %% > Auxiliary functions.
        % >> Set face f, df and d2f.
        function [fn] = Set_fn(msh)
            % >> Symbolic variable (x).
            syms x;
            % >> Deal local variables.
            [WB,EB,V,G,f_W,f_E,df_W,df_E] = ...
                deal(obj.West_Boundary,obj.East_Boundary,obj.V,obj.Gamma,obj.Phi_W,obj.Phi_E,obj.gradPhi_W,obj.gradPhi_E); 
            [x_W,x_E] = deal(msh.Xv(1),msh.Xv(msh.NV));
            L         = x_E-x_W;
            
            %  > Phi.
            if strcmpi(obj.f_Select,'Analytic')
                %  > W/o source term.
                PeL  = SubClass_2_1.Compute_Peclet;
                if     strcmpi(WB,'Dirichlet') && strcmpi(EB,'Dirichlet')
                    %  > WB: Dirichlet.
                    %  > EB: Dirichlet.
                    fn.f = f_W+((f_E-f_W).*exp(PeL.*(x./L-x_W./L))-1)./(exp(PeL)-1);
                elseif strcmpi(WB,'Dirichlet') && strcmpi(EB,'Neumann')
                    %  > WB: Dirichlet.
                    %  > EB: Neumann.
                    fn.f = f_W+df_E.*(exp(V./G.*(x-x_E))-exp(-PeL));
                elseif strcmpi(WB,'Neumann')   && strcmpi(EB,'Dirichlet')
                    %  > WB: Neumann.
                    %  > EB: Dirichlet.
                    fn.f = f_E+df_W.*(exp(V./G.*(x-x_W))-exp(PeL));
                elseif strcmpi(WB,'Dirichlet') && strcmpi(EB,'Robin')
                    %  > WB: Dirichlet.
                    %  > EB: Robin.
                    fn.f = f_W+(G./V.*df_E-f_E).*(exp(V./G.*(x-x_E))-exp(-PeL));
                elseif strcmpi(WB,'Robin')     && strcmpi(EB,'Dirichlet')
                    %  > WB: Robin.
                    %  > EB: Dirichlet.
                    fn.f = f_E+(G./V.*df_W-f_W).*(exp(V./G.*(x-x_W))-exp(PeL));
                end
            elseif strcmpi(obj.f_Select,'Manufactured')
                %  > W/ source term.
                sig  = 0.10;
                mu   = 0.50;
                fn.f = 1./(sig.*sqrt(2.*pi)).*exp(-1./2.*((x-mu)./sig).^2);
            end
            %  > gradPhi.
            fn.df   = diff(fn.f,x);
            %  > lapPhi.
            fn.d2f  = diff(diff(fn.f,x),x);
            %  > func.
            fn.func = V.*fn.df-G.*fn.d2f;
        end
        
        % >> Evaluate at...
        function [bnd,blk] = Compute_f_df_d2f(fn,msh)
            % >> Symbolic variable(x).
            syms x;
            
            % >> Evaluate f and df...
            %  > Boundaries.
            for i = 1:msh.NC+1
                bnd.f  (i) = double(subs(fn.f  ,x,msh.Xv(i)));
                bnd.df (i) = double(subs(fn.df ,x,msh.Xv(i)));
                bnd.d2f(i) = double(subs(fn.d2f,x,msh.Xv(i)));
            end
            %  > Bulk.
            for i = 1:msh.NC
                blk.f  (i) = double(subs(fn.f  ,x,msh.Xc(i)));
                blk.df (i) = double(subs(fn.df ,x,msh.Xc(i)));
                blk.d2f(i) = double(subs(fn.d2f,x,msh.Xc(i)));
            end
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
                I_i (i) = wg(i).*J.*double(subs(f,x,F_xi(i)));
            end
            I = sum(I_i);
        end

        % >> #3.
        function [F_Vol] = Compute_SourceTerm(n,fn,msh)
            for i = 1:msh.NC
                %  > Exact value (NOT USED).
                F_Vol.Ex(i) = SubClass_2_1.ExactIntegral(fn.func,msh.Xv(i),msh.Xv(i+1));
                %  > Approximated value.
                F_Vol.Ap(i) = SubClass_2_1.ApproxIntegral(n,fn.func,msh.Xv(i),msh.Xv(i+1));
            end
            F_Vol.Ex = sparse(F_Vol.Ex);
            F_Vol.Ap = sparse(F_Vol.Ap);
        end
    end
end