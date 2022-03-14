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
        %  > Estimate higher-order derivatives.
        function [pde] = T1(inp,msh,pde,s)
            %  > Stencil(s)/PDE solution(s).
            [m,~] = size(inp.ee.p);
            stl   = B_2_3_1D.SetUp_stl(inp,msh);
            for i = 1:m
                sn{i} = B_2_1_1D.Update_s(inp,msh,pde,s,stl{i});
                if i ~= m
                    xc(:,i)     = B_2_1_1D.Update_xc(sn{i});
                else
                    tt{i}       = linspace(sn{i}.p(1,i)+1,sn{i}.p(2,i),diff(sn{i}.p(:,i)));
                    sn{i}.nt(i) = length(tt{i});
                end
            end
            %  > Analytic/PDE approximated derivatives.
            v        = B_2_1_1D.Update_xv  (sn{m},xc(:,m-1)); 
            pde.df.a = B_2_1_1D.Compute_dfA(sn{m},msh.f.Xv,pde.fn.f{1});
            pde.df.x = B_2_1_1D.Compute_dfN(sn{m},tt,v); 
            %  > Plot...
        end 
        
        %% > 3. -----------------------------------------------------------
        %  > Estimate truncation error.
        function [pde] = T2(inp,msh,pde,s)
            % >> Stencil(s)/PDE solution(s).
            [m,~] = size(inp.ee.p);
            stl   = B_2_3_1D.SetUp_stl(inp,msh);
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
            %  > Update 'pde' structure.
            pde.et.av = et_a;
            pde.et.df = df;
            pde.et.fv = fv;           
            %  > Plot...
            Fig_4_1D.Plot(inp,msh,pde);
        end
    end
end