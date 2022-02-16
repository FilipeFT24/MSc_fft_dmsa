classdef B_2_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Select faces for convective/diffusive p-refinement based on exact/estimated local truncation error.
        %  > Selection criterion: local mean error (norm-1).
        function [stl] = Select_f_CD(ao,stl,ecc,et)
            [m,n] = size(et.f);
            for i = 1:n
                s{i} = find(et.f(:,i) > et.n.c(1));
            end
            stl = B_2_2_1D.Increase_p(ao,m,n,stl,s,ecc);
            % >> Rule check.
            %  > #1.
            stl = B_2_2_1D.Rule_1(ao,m,n,stl,ecc);
            %  > #2.
            
        end
        %  > 1.1.1. -------------------------------------------------------
        %  > Increase method's order.
        function [stl] = Increase_p(ao,m,n,stl,s,ecc)
            for i = 1:n
                stl.s{i} = s{i};
                if isempty(s{i})
                    continue;
                else
                    for j = 1:m
                        if ismembc(j,s{i})
                            %  > Check and increase...
                            [stl.p{i}(j,1),stl.t{i}(j,1)] = ...
                                B_2_2_1D.Check_t(ao,stl.p{i}(j,1),stl.t{i}(j,1),j,m,ecc);
                        end
                    end
                end
            end
        end
        %  > 1.1.2. -------------------------------------------------------
        %  > Check method's order and increase it accordingly.
        function [p,t] = Check_t(ao,stl_p,stl_t,j,m,ecc)
            %  > Allow uncentered differencing schemes(?).
            switch ao
                case false
                    %  > Allow CDS only.
                    switch stl_t
                        case "CDS"
                            p = stl_p+1;
                        otherwise
                            p = stl_p;
                    end
                    t = "CDS";
                case true
                    %  > Allow UDS/CDS/DDS.
                    switch stl_t
                        case "CDS"
                            p = stl_p;
                            if j-p < 1
                                t = "DDS";
                            elseif j+p+1 > m-1
                                t = "UDS";
                            else
                                if ecc(j-p) < ecc(j+p+1)
                                    t = "UDS";
                                else
                                    t = "DDS";
                                end
                            end
                            p = p+1;
                        otherwise
                            p = stl_p;
                            t = "CDS";
                    end
                otherwise
                    return;
            end
        end      
        %  > 1.1.3. -------------------------------------------------------
        %  > Compute method's order...
        function [p] = Compute_p(stl_p,stl_t)
            switch stl_t
                case "CDS"
                    p = 2.*stl_p;
                otherwise
                    p = 2.*stl_p-1;
            end
        end
        %  > 1.1.4. -------------------------------------------------------
        %  > Check rule #1...
        function [stl] = Rule_1(ao,m,n,stl,ecc)  
            for i = 1:n
                %  > Compute method's order.
                for j = 1:m
                    p(j,i) = B_2_2_1D.Compute_p(stl.p{i}(j),stl.t{i}(j));
                end
                %  > Check and increase...
                for j = 1:m
                    if j == 1 || j == m
                        continue;
                    else
                        if (p(j,i)+1 <= p(j-1,i) && p(j,i)+1 <= p(j+1,i)) %|| (p(j,i)+1 < p(j-1,i)) || (p(j,i)+1 < p(j+1,i)) 
                            [stl.p{i}(j,1),stl.t{i}(j,1)] = ...
                                B_2_2_1D.Check_t(ao,stl.p{i}(j,1),stl.t{i}(j,1),j,m,ecc);
                            stl.s{i} = cat(1,stl.s{i},j);
                            p  (j,i) = B_2_2_1D.Compute_p(stl.p{i}(j,1),stl.t{i}(j,1));
                        end
                    end
                end
                stl.s{i} = sort(stl.s{i});
            end
        end
    end
end