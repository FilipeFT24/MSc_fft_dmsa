classdef B_2_2_1D
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        %  > Select faces for convective/diffusive p-refinement based on exact/estimated local face truncation error.
        %  > Selection criterion: local mean face truncation error.
        function [stl] = Adapt_p(obj,stl,e)
            %  > Auxiliary variables.
            [m,n] = size(e.t.f);
            
            %  > Select faces for coarsening(fc)/refinement(fr) and decrease/increase method's order accordingly.
            for i = 1:n-1
                fr{i} = find(e.t.f_abs(:,i) > e.t.n_abs.f(n));
            end
            stl = B_2_2_1D.Change_p(obj,stl,fr);
            %  > Group faces and update 'stl' structure.
            if obj.gf
                %  > Method's order.
                for i = 1:n-1
                    for j = 1:m
                        p(j,i) = A_2_1D.Compute_p(stl.p{i}(j),stl.t{i}(j));
                    end
                end
                %  > Update 'stl' structure.
                stl = B_2_2_1D.Group_f(obj,stl,p);
            end
        end
        % >> 1.2. ---------------------------------------------------------
        %  > Decrease/increase method's order accordingly.
        function [stl] = Change_p(obj,stl,fr)
            for i = 1:size(fr,2)
                if isempty(fr{i})
                    stl.s{i} = [];
                else
                    l = 0;
                    for j = 1:size(fr{i},1)
                        %  > Auxiliary variables.
                        k = fr{i}(j);
                        l = l+1;
                        
                        %  > Increase method's order.
                        [stl.p{i}(k),stl.t{i}(k)] = ...
                            B_2_2_1D.Increase_p(obj,stl.p{i}(k),stl.t{i}(k));
                        %  > Add element to 'stl.s'.
                        stl_s{i}(l,1) = k;
                    end
                    stl.s{i} = stl_s{i};
                end
            end
        end
        %  > 1.2.1. -------------------------------------------------------
        %  > Increase method's order.
        function [stl_p,stl_t] = Increase_p(obj,stl_p,stl_t)
            %  > Allow uncentered differencing schemes(?).
            switch obj.ao
                case false
                    %  > Allow CDS only.
                    switch stl_t
                        case "CDS"
                            stl_p = stl_p+1;
                        otherwise
                            stl_t = "CDS";
                    end
                case true
                    %  > Allow UDS/CDS/DDS.
                    switch stl_t
                        case "CDS"
                            stl_p = stl_p+1;
                            stl_t = "UDS";
                        otherwise
                            stl_t = "CDS";
                    end
            end
        end
        %  > 1.2.2. -------------------------------------------------------
        %  > Decrease method's order.
        function [] = Decrease_p()
        end
        % >> 1.3. ---------------------------------------------------------
        %  > Group faces.
        function [stl] = Group_f(obj,stl,p)
            %  > Auxiliary variables.
            [m,n] = size(p);
            
            
%             hold on;
%             plot(p(:,1),'-m');
%             plot(p(:,2),'-c');
%             
            
            %  > Group by order.
            for i = 1:n
                k = 0;
                for j = 1:m
                    if (j > 1) && (p(j,i) == p(j-1,i))
                        l            = l+1;
                        g{i}{k}(l,1) = j;
                    else
                        k            = k+1;
                        l            = 1;
                        g{i}{k}(l,1) = j;
                    end
                end
            end
            %  > Increase method's order.
            for i = 1:n
                j = size(g{i},2);
                if j ~= 1
                    for k = 1:j
                        v{i}(k,1) = p(g{i}{k}(1),i);
                    end
                else
                    v{i} = p(g{i}{j}(1),i);
                end
            end
            %  > Select refinement groups.
            for i = 1:n
                j = size(g{i},2);
                if j ~= 1
                    for k = 1:j
                        if k == 1 || k == j
                            continue;
                        else
                            if (v{i}(k) < v{i}(k-1)) && (v{i}(k) < v{i}(k+1))
                                for l = 1:size(g{i}{k},1)
                                    %  > Auxiliary variables.
                                    s = size(stl.s{i},1);
                                    o = g{i}{k}(l);
                                    
                                    %  > Increase method's order (update fields).
                                    [stl.p{i}(o),stl.t{i}(o)] = ...
                                        B_2_2_1D.Increase_p(obj,stl.p{i}(o),stl.t{i}(o));
                                    stl.s{i}(s+1) = o;
                                    %  > Update vector p (redundant).
                                    p(o,i) = A_2_1D.Compute_p(stl.p{i}(o),stl.t{i}(o));
                                end
                            end
                        end
                    end
                end
                stl.s{i} = sort(stl.s{i});
            end
            
%             hold on;
%             plot(p(:,1),'-r');
%             plot(p(:,2),'-b');
            
            
        end
    end
end