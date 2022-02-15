classdef B_2_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Select faces for convective/diffusive p-refinement based on local truncation error...
        function [stl] = Select_f_CD(et,stl,allow_odd)
            % >> Increase method's order based on mean error (norm-1).
            %  > Size.
            [m,n] = size(et.f);
            for i = 1:n
                s{i} = find(et.f(:,i) > et.n.c(1));
            end
            stl = B_2_2_1D.Increase_p(s,stl,allow_odd,m,n);
            % >> Rule check...
            stl = B_2_2_1D.Rule_1(stl,allow_odd,m,n);
        end
        %  > 1.1.1. -------------------------------------------------------
        %  > Increase method's order...
        function [stl] = Increase_p(s,stl_pst,allow_odd,m,n)
            %  > Allow uncentered differencing schemes(?).
            for i = 1:n
                stl.s{i} = s{i};
                for j = 1:m
                    if ~ismembc(j,s{i})
                        stl.p{i}(j,1) = stl_pst.p{i}(j);
                        stl.t{i}(j,1) = stl_pst.t{i}(j);
                    else
                        switch allow_odd
                            case false
                                %  > Allow CDS only...
                                stl.p{i}(j,1) = stl_pst.p{i}(j)+1;
                                stl.t{i}(j,1) = stl_pst.t{i}(j);
                            case true
                                %  > Allow UDS/CDS/DDS...
                            otherwise
                                return;
                        end
                    end
                end
            end
        end
        %  > 1.1.2. -------------------------------------------------------
        %  > Compute method's order...
        function [p] = Compute_p(stl_p,stl_t)
            switch stl_t
                case "CDS"
                    p = 2.*stl_p;
                otherwise
                    p = 2.*stl_p-1;
            end
        end
        %  > 1.1.3. -------------------------------------------------------
        %  > Check rule #1...
        function [stl] = Rule_1(stl,allow_odd,m,n)
            for i = 1:n
                for j = 1:m
                    p(j,i) = B_2_2_1D.Compute_p(stl.p{i}(j),stl.t{i}(j));
                end
                %  > Increase...
                for j = 1:m
                    if j == 1 || j == m
                        continue;
                    else
                        if p(j,i)+1 < p(j-1,i) && p(j,i)+1 < p(j+1,i)
                            switch allow_odd
                                case false
                                    %  > Allow CDS only...
                                    stl.s{i}      = cat(1,stl.s{i},j);
                                    stl.p{i}(j,1) = stl.p{i}(j)+1;
                                    p       (j,i) = B_2_2_1D.Compute_p(stl.p{i}(j,1),stl.t{i}(j,1));
                                case true
                                    %  > Allow UDS/CDS/DDS...
                                otherwise
                                    return;
                            end
                        end
                    end
                end
                stl.s{i} = sort(stl.s{i});
            end
        end
    end
end