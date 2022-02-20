classdef B_2_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Select faces for convective/diffusive p-refinement based on exact/estimated local truncation error.
        %  > Selection criterion: local mean error (norm-1).
        function [stl] = Select_fCD(ao,stl,ecc,et)
            [m,n] = size(et.f);
            %  > Step #1: Increase method's order.
            for i = 1:n
                s{i} = find(et.f(:,i) > et.n.t(1));
            end
            stl = B_2_2_1D.Increase_p(ao,m,n,stl,s,ecc);
            %  > Step #2: Check rules.
            R1  = true;
            R2  = true;
            stl = B_2_2_1D.Check_Rules(R1,R2,ao,m,n,stl,ecc);
        end
        %  > 1.1.1. -------------------------------------------------------
        %  > Increase method's order.
        function [stl] = Increase_p(ao,m,n,stl,s,ecc)
            for i = 1:n
                stl.s{i} = s{i};
                if ~isempty(s{i})
                    for j = 1:m
                        if ismembc(j,s{i})
                            %  > Update stl.
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
        %  > Compute method's order.
        function [p] = Compute_p(stl_p,stl_t)
            switch stl_t
                case "CDS"
                    p = 2.*stl_p;
                otherwise
                    p = 2.*stl_p-1;
            end
        end
        %  > 1.1.4. -------------------------------------------------------
        %  > Group cell faces (...to apply rule #1).
        function [pc] = Group_CellFaces(stl,m,n)
            %  > Compute method's order.
            for i = 1:n
                for j = 1:m
                    pf(j,i) = B_2_2_1D.Compute_p(stl.p{i}(j),stl.t{i}(j));
                end
            end
            %  > Group...
            for i = 1:m-1
                for j = 1:n
                    pc{i,1}(j,:) = [pf(i,j),pf(i+1,j)];
                end
            end
        end
        %  > 1.1.5. -------------------------------------------------------
        %  > Check rule #1...
        function [stl] = Rule_1(ao,m,n,stl,ecc,pc)
            for i = 1:n
                for j = 1:m-1
                    if j == 1 || j == m-1
                        continue;
                    else
                        %  > Evaluated cell (j).
                        je = pc{j,1}(i,:) ;
                        %  > Minimum elements of cells (j-1) and (j+1).
                        jm = max(pc{j-1,1}(i,:));
                        jp = max(pc{j+1,1}(i,:));
                        %  > Check...
                        if all(je < jm) && all(je < jp)
                            %  > Update stl.
                            [stl.p{i}(j),stl.t{i}(j)] = ...
                                B_2_2_1D.Check_t(ao,stl.p{i}(j,1),stl.t{i}(j,1),j,m,ecc);
                            stl.s{i} = cat(1,stl.s{i},j);
                            %  > Update pc (redundant).
                            pc{j-1,1}(i,2) = stl.p{i}(j);
                            pc{j  ,1}(i,1) = stl.p{i}(j);
                        end
                    end
                end
                stl.s{i} = sort(stl.s{i});
            end
        end   
        %  > 1.1.6. -------------------------------------------------------
        %  > Check rule #2...
        function [stl] = Rule_2(ao,m,n,stl,ecc)  
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
                        if p(j,i)+1 <= p(j-1,i) && p(j,i)+1 <= p(j+1,i)
                            %  > Update stl.
                            [stl.p{i}(j,1),stl.t{i}(j,1)] = ...
                                B_2_2_1D.Check_t(ao,stl.p{i}(j,1),stl.t{i}(j,1),j,m,ecc);
                            stl.s{i} = cat(1,stl.s{i},j);
                            %  > Update p.
                            p  (j,i) = B_2_2_1D.Compute_p(stl.p{i}(j,1),stl.t{i}(j,1));
                        end
                    end
                end
                stl.s{i} = sort(stl.s{i});
            end
        end
        %  > 1.1.7. -------------------------------------------------------
        %  > Check rules...
        function [stl] = Check_Rules(R1,R2,ao,m,n,stl,ecc)
            %  > Rule #1.
            switch R1
                case true
                    %  > Group...
                    pc  = B_2_2_1D.Group_CellFaces(stl,m,n);
                    stl = B_2_2_1D.Rule_1(ao,m,n,stl,ecc,pc);
                otherwise
                    %  > Do not check rule #1...
            end
            %  > Rule #2.
            switch R2
                case true
                    stl = B_2_2_1D.Rule_2(ao,m,n,stl,ecc);
                otherwise
                    %  > Do not check rule #2...
            end
        end
    end
end