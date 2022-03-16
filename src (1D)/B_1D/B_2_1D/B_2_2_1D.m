classdef B_2_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        %  > Select faces for convective-diffusive p-refinement.
        function [pde,stl] = Adapt_p(inp,pde,stl)
            [~,n] = size(pde.e.t.p);
            for i = 1:n
                m(i) = 1;
                f{i} = find(pde.e.t.p{m(i)}.f_abs(:,i) > pde.e.t.p{m(i)}.n_abs.f(n));
            end
            stl = B_2_2_1D.Update_stl(inp,stl,f);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Increase method's order.
        function [stl] = Update_stl(inp,stl,f)
            for i = 1:size(f,2)
                n = 2;
                if ~isempty(f{i})
                    for j = 1:size(f{i},1)
                        k = i*size(f,2)-1;
                        if ~inp.pa.odd
                            l                = n;
                            stl.p(f{i}(j),k) = inc(stl.p(f{i}(j),k),l).i;
                        end
                    end
                end
                stl.s{i} = f{i};
            end
        end
    end
end