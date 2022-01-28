classdef B_Tools
    methods (Static)
        %% > B_Tools.
        % >> --------------------------------------------------------------
        % >> --------------------------------------------------------------
        
        %% > 1. -----------------------------------------------------------
        function Y = LU(X,str)
            %  > Initialization.
            [m,n]     = size(X);
            [U,L,B,Y] = deal(zeros(n,n));
            
            if strcmpi(str,'Dolittle')
                L = ones(n,n);
                for i = 1:n
                    j      = i:n;
                    sum    = 0;
                    k      = 1:i-1;
                    sum    = sum + L(i,k)*U(k,j);
                    U(i,j) = X(i,j)-sum;
                    sum    = 0;
                    k      = 1:i-1;
                    sum    = sum + L(j,k)*U(k,i);
                    L(j,i) = (X(j,i)-sum)/U(i,i);
                end
            elseif strcmpi(str,'Crout')
                for j = 1:n
                    for i = 1:n
                        sum = 0;
                        k   = 1:j-1;
                        sum = sum+(L(i,k)*U(k,j));
                        if i == j
                            U(i,j) = 1;
                        end
                        if i >= j
                            L(i,j) = X(i,j)-sum;
                        else
                            U(i,j) = (1/L(i,i))*(X(i,j)-sum);
                        end
                    end
                end
            end
            %  > Forward elimination (LB=I).
            b      = eye(n);
            i      = 1:m;
            B(1,i) = b(1,i)/L(1,1);
            for k  = 2:m
                sum    = 0;
                j      = k-1:-1:1;
                sum    = sum+L(k,j)*B(j,i);
                B(k,i) = (b(k,i)-sum)/L(k,k);
            end
            %  > Backward substitution (U*A^-1 = B).
            i      = 1:m;
            Y(m,i) = B(m,i)/U(m,m);
            for k  = m-1:-1:1
                sum    = 0;
                j      = k+1:m;
                sum    = sum+U(k,j)*Y(j,i);
                Y(k,i) = (B(k,i)-sum)/U(k,k);
            end
        end
    end
end