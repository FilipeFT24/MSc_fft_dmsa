classdef SubClass_1_3_3
    methods (Static)
        %% > SubClass_1_3_3.
        % >> --------------------------------------------------------------
        % >> 1.   Tools.
        %  > 1.1. Modified version(#1) of ismember/ismembc: Compare array  B w/ sorted array  A.
        %  > 1.2. Modified version(#2) of ismember/ismembc: Compare array  B w/ sorted array  A and return nnz elements of B if nnz(Array) > 2.
        %  > 1.3. Modified version(#1) of isequal: Compare array   B w/ matrix A.
        %  > 1.4. Modified version(#2) of isequal: Compare matrix  B w/ matrix A.
        %  > 1.5. Modified version(#1) of setdiff.
        % >> --------------------------------------------------------------

        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [Flag] = fft_ismember_1(A,B)
            Flag = double(any(ismembc(B,sort(A))));
        end
        % >> 1.2. ---------------------------------------------------------
        function [Flag] = fft_ismember_2(A,B)
            Array = find(double(ismembc(B,sort(A))));
            if size(Array) < 2
                Flag = zeros(1,2);
            else
                Flag = B(Array);
            end  
        end
        % >> 1.3. ---------------------------------------------------------
        function [i_Flag] = fft_isequal_1(A,B)
            i_Flag = 0;
            for i = 1:size(A,1)
                for j = 1:size(B,1)
                    Flag(i,j) = isequal(A(i,:),B(j,:));
                end
                if any(Flag(i,:))
                    i_Flag = 1;
                    break;
                end
            end
        end
        % >> 1.4. ---------------------------------------------------------
        function [i_Flag] = fft_isequal_2(A,B)
            i_Flag = 0;
            for i = 1:size(A,1)
                for j = 1:size(B,1)
                    Flag(i,j) = isequal(A(i,:),B(j,:));
                    continue;
                end
            end
            if nnz(Flag) > 1
                i_Flag = 1;
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
    end
end