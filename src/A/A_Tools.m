classdef A_Tools
    methods (Static)
        %% > A_Tools.
        % >> --------------------------------------------------------------
        % >> 1.   Tools.
        %  > 1.1. Modified version (#1) of ismember/ismembc: Compare array  B w/ sorted array  A.
        %  > 1.2. Modified version (#2) of ismember/ismembc: Compare array  B w/ sorted array  A and return nnz elements of B if nnz(Array) > 2.
        %  > 1.3. Modified version (#1) of isequal: Compare array   B w/ matrix A.
        %  > 1.4. Modified version (#2) of isequal: Compare matrix  B w/ matrix A.
        %  > 1.5. Modified version (#1) of setdiff.
        %  > 1.6. Modified version (#1) of 'unique'.
        %  > 1.7. Modified version (#1) of 'mean'.
        %  > 1.8. Modified version (#1) of 'pdist': Point/Point. 
        %  > 1.9. Modified version (#1) of 'pdist': Point/Matrix.
        % >> 2.   Sort 'msh' fields.
        % >> --------------------------------------------------------------
        
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [Flag] = fft_ismember_1(A,B)
            Flag = any(ismembc(B,sort(A)));
        end
        % >> 1.2. ---------------------------------------------------------
        function [Flag] = fft_ismember_2(A,B)
            Array = find(double(ismembc(B,A)));
            if size(Array) < 2
                Flag = zeros(1,2);
            else
                Flag = B(Array);
            end
        end
        % >> 1.3. ---------------------------------------------------------
        function [i_Flag] = fft_isequal_1(A,B)
            i_Flag = false;
            for i = 1:size(A,1)
                for j = 1:size(B,1)
                    Flag(i,j) = isequal(A(i,:),B(j,:));
                end
                if any(Flag(i,:))
                    i_Flag = true;
                    break;
                end
            end
        end
        % >> 1.4. ---------------------------------------------------------
        function [i_Flag] = fft_isequal_2(A,B)
            i_Flag = false;
            for i = 1:size(A,1)
                for j = 1:size(B,1)
                    Flag(i,j) = isequal(A(i,:),B(j,:));
                    continue;
                end
            end
            if nnz(Flag) > 1
                i_Flag = true;
            end
        end
        % >> 1.5. ---------------------------------------------------------
        function [Z] = fft_setdiff(X,Y)
            if ~isempty(X) && ~isempty(Y)
                check    = false(1, max(max(X), max(Y)));
                check(X) = true;
                check(Y) = false;
                Z        = X(check(X));
            else
                Z = X;
            end
        end
        % >> 1.6. ---------------------------------------------------------
        function [Y] = fft_unique(X)
            Y = [X(not(~diff(X))),X(length(X))];
        end
        % >> 1.7. ---------------------------------------------------------
        function [Y] = fft_mean(X)
            Y = sum(X)./length(X);
        end
        % >> 1.8. ---------------------------------------------------------
        function [D] = fft_dist_1(A,B)
            D = sqrt(abs(bsxfun(@plus,sum(A.*A,1),sum(B.*B,1).')-2*A'*B));
        end
        % >> 1.9. ---------------------------------------------------------
        function [D] = fft_dist_2(A,B)
            D = sqrt(abs(bsxfun(@plus,sum(B.*B,1),sum(A.*A,1).')'-2*B'*A));
        end
        
        %% > 2. -----------------------------------------------------------      
        % >> 2.2. ---------------------------------------------------------
        function [msh] = Sort_msh(msh)
            % >> msh.
            msh     = orderfields(msh        ,{'d','c','f','bnd','s'});
            %  > c.
            msh.c   = orderfields(msh.c      ,{'NC','xy_v','mean','h','c','f'});
            msh.c.f = orderfields(msh.c.f    ,{'f','xy_v','mean','len','Nf','Sf'});
            %  > f.
            msh.f   = orderfields(msh.f      ,{'NF','xy_v','mean','c'});
            %  > bnd.
            msh.bnd = orderfields(msh.bnd    ,{'c','f'});
            %  > s.
            msh.s     = orderfields(msh.s    ,{'c','f','c_e','f_e','xy_v_c','xy_v_f','xy_v_t','par'});
            msh.s.par = orderfields(msh.s.par,{'ne','nx','ny','lx','ly'});
        end
    end
end