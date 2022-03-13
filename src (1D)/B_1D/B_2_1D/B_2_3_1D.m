classdef B_2_3_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        %  > Initialize 'stl' structure.
        function [stl] = SetUp_stl(inp,msh)
            [m,n] = size(inp.ee.p);
            for i = 1:m
                for j = 1:n
                    stl{i}.p   (:,j) = repelem(inp.ee.p(i,j),msh.f.NF);
                    stl{i}.s{j}(:,1) = 1:msh.f.NF;
                    stl{i}.t   (:,j) = repelem(inp.ee.s(i,j),msh.f.NF);
                end
            end
        end

        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Compute derivatives up to order n (w/ analytic solution).
        function [df] = Compute_dfa(msh,f,p,nt)
            syms x;
            for i = 1:length(p)
                for j = 0:nt(i)-1
                    n      = p(i)+j;
                    dfn{i} = diff(f,x,n);
                    dfn{i} = matlabFunction(dfn{i});
                    if nargin(dfn{i}) ~= 0
                        df{i}(:,j+1) = dfn{i}(msh.f.Xv)./factorial(n);
                    else
                        df{i}(:,j+1) = zeros(msh.f.NF,1);
                    end
                end
            end
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Compute magnitude of truncated terms (w/ analytic solution).
        function [msh,pde] = T_1(inp,msh,pde,s,nt)
            %  > Stencil
            stl     = B_2_3_1D.SetUp_stl  (inp,msh);
            s       = B_2_1_1D.Update_s   (inp,msh,pde,s,stl{:});
            msh     = Tools_1D.Set_msh    (msh,stl,s);
            %  > Truncated terms.            
            pde.e.t = B_2_1_1D.Update_et_f(inp,s,pde.av);
            p       = A_2_1D.Compute_p    (inp.ee.p,inp.ee.s);
            df      = B_2_3_1D.Compute_dfa(msh,pde.fn.f{1},p,nt);
            pde.e.m = A_2_1D.Compute_tm   (inp,msh,s,stl{:},df,p,nt);
            %  > Plot...
            Fig_2_1D.Plot(msh,pde,p);
        end

        %% > 3. -----------------------------------------------------------
        % >> 3.1. ---------------------------------------------------------
        %  > Compute derivatives w/ higher-order solution up to order n.
        function [dfn] = Compute_dfn(s,stl,x_c)
            [m,n] = size(s.Inv);
            [v,~] = B_2_1_1D.Update_xf(s,x_c);
            for i = 1:m
                j = stl.o(i)+1;
                k = stl.o(i)+stl.nt(i);
                for l = 1:n
                    for o = j:k
                        p           = o+1-j;
                        df          = zeros(1,size(s.Inv{i,l},1));
                        df    (1,o) = 1;
                        dfn{i}(l,p) = df*s.Inv{i,l}*v{i,l};
                    end
                end
            end
        end
        % >> 3.2. ---------------------------------------------------------
        %  > Truncated terms' magnitude (w/ higher-order solution).
        function [] = EE_2(obj,inp_lo,inp_ho,msh,pde,s)
            %  > Initialize.
            [stl_lo]    = B_2_3_1D.SetUp_stl  (inp_lo,msh);
            [stl_lo.nt] = inp_lo.nt;
            [stl_ho]    = B_2_3_1D.SetUp_stl  (inp_ho,msh);
            %  > 'p-standard' run.
            [s_lo]      = B_2_1_1D.Update_s   (obj,msh,pde,s,stl_lo);
            [s_ho]      = B_2_1_1D.Update_s   (obj,msh,pde,s,stl_ho);
            [x_ho.c]    = B_2_1_1D.Update_xc  (s_ho);
            %  > Derivatives' magnitude.
            [df.a]      = B_2_3_1D.Compute_dfa(msh,pde,stl_lo);
            [df.n]      = B_2_3_1D.Compute_dfn(s_ho,stl_lo,x_ho.c);
            %  > Compute truncated terms.
            [tm.a]      = A_2_1D.Compute_tm   (obj,msh,s_lo,stl_lo,df.a);
            [tm.n]      = A_2_1D.Compute_tm   (obj,msh,s_lo,stl_lo,df.n);
            %  > Plot.
            Fig_2_1D.WrapUp_Fig_2_2_1D(msh,stl_lo.o-1,df,tm); 
        end 
        
        %% > 4. -----------------------------------------------------------
        % >> 4.1. ---------------------------------------------------------
        function [pde,p] = T_3(inp,msh,pde,s)
            %  > Stencil(s)/PDE solution(s).
            m   = size(inp.ee.p);
            stl = B_2_3_1D.SetUp_stl(inp,msh);
            for i = 1:m
                sn  {i} = B_2_1_1D.Update_s (inp,msh,pde,s,stl{i});
                xc(:,i) = B_2_1_1D.Update_xc(sn{i}); 
            end
 
            % >> Face values.
            for i = 1:m
                %  > w/ analytic solution.
                fv.a{i} = B_2_1_1D.Update_xf(sn{i},pde.av.c);
                %  > w/ PDE solution: i-th order PDE solution/j-th order stencil coefficients.
                for j = 1:m
                    fv.x{i,j} = B_2_1_1D.Update_xf(sn{j},xc(:,i));
                end
            end
            %  > Approximated analytic-PDE face value. 
            for i = 1:m
                for j = 1:m
                    df.ax{i,j} = abs(fv.a{j}-fv.x{i,j});
                end
            end
            % >> Truncation error/truncation error difference.
            for i = 1:m
                %  > w/ analytic solution.
                et_a{i}   = abs(pde.av.f-fv.a{i});
                diff(i,:) = setdiff(1:m,i);
                for j = 1:m-1
                    df.a{i,j} = abs(fv.a{i}-fv.a{diff(i,j)});
                end
                %  > w/ PDE solution.
                for j = 1:m-1
                    df.x{i,j} = abs(fv.x{i,j}-fv.x{i,j+1});
                end
            end
            % >> Curve translation difference.
            for i = 1:m
                for j = 1:m-1
                    df.t{i,j} = abs((fv.a{j+1}-fv.x{i,j+1})-(fv.a{j}-fv.x{i,j}));
                end
            end
            pde.et.av = et_a;
            pde.et.df = df;
            pde.et.fv = fv;
            p         = A_2_1D.Compute_p(inp.ee.p,inp.ee.s);
            %  > Plot...
            Fig_4_1D.Plot(msh,pde,p);
        end
    end
end