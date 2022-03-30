classdef B_2_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > 'p-standard' run.
        function [obj,msh] = p_standard(inp,obj,msh,pde)
            % >> Update fields 'm', 's' and 'x'.
            add          = 1;
            obj.s        = B_2_1_1D.Update_1   (inp,msh,pde,obj.s,obj.u,add);
            obj.x        = B_2_1_1D.Update_2   (inp,msh,obj.s,obj.u,obj.x);
            obj.m        = B_2_1_1D.Update_3   (msh,pde,obj.m,obj.s,obj.u,obj.x);
            % >> Update...
            %  > ...solution.
            obj.x.nv.a   = pde.av;
            obj.x.nv.x.c = B_2_1_1D.Update_xc  (obj.m.At,obj.m.Bt);
            %  > ...field 'x'.
            obj.x        = B_2_1_1D.Update_4   (obj.s,obj.u,obj.x);
            % >> Update cell/face truncation and cell global discretization error distribution/norms.
            obj.e        = B_2_1_1D.Update_e   (inp,msh,pde,obj.e,obj.m,obj.s,obj.u,obj.x,0);
            % >> Compute tuncated terms' magnitude(?).
            if inp.pt.tt
                obj.e.t = B_2_2_1D.p_truncation(inp,msh,pde,obj.e.t,obj.s,obj.u,obj.x);
            end
            % >> Update structures.
            [obj,msh]    = B_2_1_1D.Set_struct (obj,msh);
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Check truncated terms' magnitude.
        function [et] = p_truncation(inp,msh,pde,et,s,u,x)
            %  > Compute analytic derivatives (dfA).
            dfA = B_2_2_1D.Compute_dfA(inp,pde.fn.f{1});
            %  > Compute truncated terms' magnitude.
            et  = B_2_2_1D.Compute_TTM(inp,msh,dfA,et,s,u,x);
        end
        %  > 1.2.1. -------------------------------------------------------
        %  > Compute derivatives (w/ analytic solution).
        function [dfA] = Compute_dfA(inp,f)
            %  > Auxiliary variables.
            syms x;
            n = inp.pt.nt;
            p = inp.ps.p;
            q = length(p);
           
            %  > Derive up to...
            for i = 1:q
                j(i) = p(i)+n(i);
            end
            %  > Derivative(s)' function handle(s)'.
            for i = 1:max(j)
                dfA{i} = diff(f,x,i)/factorial(i);
                dfA{i} = matlabFunction(dfA{i});
            end
        end
        %  > 1.2.2. ---------------------------------------------------------
        %  > Assemble matrix Df' (w/ analytic solution).
        function [et] = Compute_TTM(inp,msh,dfn,et,s,u,x)
            %  > Auxiliary varibales.
            nt = inp.pt.nt;
            ls = size(u.s,2);
            
            for i = 1:ls
                if isempty(u.s{i})
                    continue;
                else
                    for j = 1:size(u.s{i},1)
                        k = u.s{i}(j);
                        X = s.t{k,i};
                        F = msh.f.Xv(k);
                        
                        %  > Df.
                        df = zeros(1,nt(i));
                        Df = zeros(length(X),nt(i));
                        switch string(s.bt{k,i})
                            case "Neumann"
                                l         = 1:nt(i);
                                m         = 1:length(X);
                                n         = 1:length(X)-1;
                                o         = stl.o(i)-1+l;
                                q         = n-1;
                                r         = length(X);
                                t         = o-1;
                                Df  (n,l) = (X(n)-F)'.^o;
                                Df  (r,l) = (X(r)-F)'.^t.*o;
                                Df_T(l,m) = transpose(Df(m,l));
                            case "Robin"
                                g_v       = inp.g./inp.v;
                                l         = 1:nt(i);
                                m         = 1:length(X);
                                n         = 1:length(X)-1;
                                o         = s.p(i)-1+l;
                                q         = n-1;
                                r         = length(X);
                                t         = o-1;
                                Df  (m,l) = (X(m)-F)'.^o;
                                Df  (r,l) = Df(r,l)+g_v.*(X(r)-F)'.^t.*o;
                                Df  (l,m) = transpose(Df(m,l));
                            otherwise
                                l           = 1:nt(i);
                                m           = 1:length(X);
                                n           = u.p(k,i*ls-1)+l-1+i-1;
                                Df  (m,l)   = (X(m)-F)'.^n;
                                for i_df    = 1:nt(i)
                                    df(1,i_df) = dfn{n(i_df)}(F);
                                end                              
                        end
                        TTM{i}(k,l) = x.cf{k,i}*Df.*df;
                    end
                end         
                et.f    {i} = TTM{i};
                et.f_abs{i} = abs(et.f{i});
            end
            a        = 1:msh.c.Nc;
            b        = a+1;
            sum_ttm  = s.v(1).*TTM{1}+s.v(2).*TTM{2};
            et.c     = sum_ttm(b,:)-sum_ttm(a,:);
            et.c_abs = abs(et.c);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > 'p-adaptative' run.
        function [obj,msh] = p_adaptative(inp,obj,msh,pde)
            i = 0;
            while 1
                % >> Update cycle count...
                i = i+1;
                % >> Update fields 'm', 's' and 'x'.
                obj.s        = B_2_1_1D.Update_1       (inp,msh,pde,obj.s,obj.u);
                obj.x        = B_2_1_1D.Update_2       (inp,msh,obj.s,obj.u,obj.x);
                obj.m        = B_2_1_1D.Update_3       (msh,pde,obj.m,obj.s,obj.u,obj.x);
                % >> Update...
                %  > ...solution(?).
                flag_1           = B_2_2_1D.Solve_AX_B (i);
                if flag_1
                    obj.x.nv.a   = pde.av;
                    obj.x.nv.x.c = B_2_1_1D.Update_xc  (obj.m.At,obj.m.Bt);
                end
                %  > ...field 'x'.
                obj.x            = B_2_1_1D.Update_4   (obj.s,obj.u,obj.x);
                % >> Update cell/face truncation and cell global discretization error distribution/norms.
                obj.e            = B_2_1_1D.Update_et_a(msh,obj.e,obj.s,obj.u,obj.x);
                obj.e            = B_2_1_1D.Update_ec_a(msh,obj.e,obj.x);
                obj.e            = B_2_1_1D.Update_et_p(inp,msh,pde,obj.e,obj.s,obj.u,obj.x);
                obj.e            = B_2_1_1D.Update_ec_p(msh,obj.e,obj.m);
                % >> Set structures.
                obj_m(i,:)       = obj.m;
                obj_e(i,:)       = obj.e.p;
                %  >> Stop adaptation(?).
                flag_2           = B_2_2_1D.Stop       (i,obj_e{i,1}.c.n_abs);
                if ~flag_2
                    obj.u        = B_2_2_1D.Update_u   (obj_e{i,1},obj.u);
                else
                    obj.m        = obj_m;
                    obj.e        = obj_e;
                    break;
                end
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > 2.2.1. -------------------------------------------------------
        %  > Select faces for p-coarsening/refinement.
        %  > nf    : Number of selected faces.
        %  > lambda: Parameter used to establish a treshold for the maximum face truncation error (group selection).        
        function [u_s] = Select_f(e)
            %  > Auxiliary variables.
            nf     = 2;
            lambda = 0.95;
            nv     = size(e.t.f_abs,2);
            
            [~,iM] = maxk(e.t.f_abs(:,1:nv-1),nf);
            for i  = 1:nv-1
                u_s{i}(:,1) = sort(iM(:,i));
            end
        end
        %  > 2.2.2. -------------------------------------------------------
        %  > Decrease/increase method's order (update field 'u').
        function [u] = Update_u(e,u)
            %  > Auxiliary variables.
            u_s = B_2_2_1D.Select_f(e);
            m   = size(u_s,2);
            A   = 2;
                       
            for i = 1:m
                if ~isempty(u_s{i})
                    for j = 1:size(u_s{i},1)
                        %  > Face/column index.
                        k        = u_s{i}(j);
                        l        = i*m-1;
                        %  > Increase polynomial order by "A".
                        u.p(k,l) = u.p(k,l)+A;
                    end
                end
                u.s{i} = u_s{i};
            end
        end
        % >> 2.3. ---------------------------------------------------------
        %  > Solve AX=B(?) criterion/criteria.
        function [flag] = Solve_AX_B(i)
            if i ~= 1
                flag = 0;
            else
                flag = 1;
            end
        end
        % >> 2.4. ---------------------------------------------------------
        %  > Set stopping criterion/criteria.
        %  > nc_M   : Maximum number of cycles.
        %  > ec_m_L1: Minimum cell global discretization error (L1 norm).
        %  > ec_M_L3: Maximum cell global discretization error (L_infinity norm).
        function [flag] = Stop(c,ec_abs)
            %  > Auxiliary variables.
            nc_M = 3;
            ec_m = 1e-10;
            
            if c > nc_M || ec_abs(1) <= ec_m
                flag = 1;
            else
                flag = 0;
            end
        end
    end
end