classdef B_2_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        %  > Select faces for convective/diffusive p-refinement based on exact/estimated local face truncation error.
        %  > Selection criterion: local mean face truncation error.
        function [stl] = Adapt_p(inp,stl,e)
            %  > Auxiliary variables.
            [~,n] = size(e.t.f);
            
            %  > Select faces for refinement(fr) and increase method's order accordingly.
            for i = 1:n-1
                fr{i} = find(e.t.f_abs(:,i) > e.t.n_abs.f(n));
            end
            stl = B_2_2_1D.Update_stl(inp,stl,fr);
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        %  > Increase method's order accordingly.
        function [stl] = Update_stl(inp,stl,fr)
            for i = 1:size(fr,2)
                if isempty(fr{i})
                    s_fr{i} = [];
                else
                    l = 0;
                    for j = 1:size(fr{i},1)
                        %  > Auxiliary variables.
                        k = fr{i}(j);
                        l = l+1;
                        %  > Increase method's order.
                        [stl.p(k,i),stl.t(k,i)] = ...
                            B_2_2_1D.Change_p(inp,stl.p(k,i),stl.t(k,i));
                        %  > Add element to 's_fr'.
                        s_fr{i}(l,1) = k;
                    end
                end
            end
            stl.s = s_fr;
        end
        % >> 2.2. ---------------------------------------------------------
        %  > Decrease/increase method's order.
        function [stl_p,stl_t] = Change_p(obj,stl_p,stl_t)
            %  > Allow uncentered differencing schemes(?).
            switch obj.pa.odd
                case false
                    %  > Allow CDS only.
                    switch stl_t
                        case "c"
                            stl_p = stl_p+1;
                        otherwise
                            stl_t = "c";   
                    end
                case true
                    %  > Allow UDS/CDS/DDS.
            end
        end
    end
end