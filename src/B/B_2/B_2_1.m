classdef B_2_1
    methods (Static)
        %% > Wrap-up B_2_1.
        function [pde] = WrapUp_B_2_1(pde,np,wf)
            % >> ----------------------------------------------------------
            % >> 1.   Compute polynomial terms for phi and grad(phi).
            %  > 1.1. Auxiliary function (based on 'iD' input).
            % >> 2.   Set weighting function.
            % >> ----------------------------------------------------------
            % >> 1.
            [pde.pr.numb,pde.pr.Coeff_1,pde.pr.Coeff_2,pde.pr.Exp_1,pde.pr.Exp_2] = ...
                B_2_1.Polynomial_Reconstruction(np);
            % >> 2.
            if strcmpi(wf,'Weighted')
                [pde.wf.wf_1,pde.wf.wf_2] = B_2_1.W_Function();
            end
        end
        
        %% > 1. -----------------------------------------------------------
        function [len_p,Coeff_1,Coeff_2,Exp_1,Exp_2] = Polynomial_Reconstruction(p)
            %  > Number of terms.
            len_p = 1./2.*(p+1).*(p+2);
            %  > Polynomial coefficients/exponents.
            [Coeff_1,Exp_1] = B_2_1.Compute_PolTerms(1,p,len_p);
            [Coeff_2,Exp_2] = B_2_1.Compute_PolTerms(2,p,len_p);
        end
        % >> 1.1. ---------------------------------------------------------
        function [Coeff_iD,Exp_iD] = Compute_PolTerms(iD,p,len_p)
            %                |------------------------------------------------------------------------------------------------|
            %  > Term      = | 1 | x | y | x^2 | xy  | y^2 | x^3  | x^2y | xy^2 | y^3  | x^4  | x^3y  | x^2y^2 | xy^3  | y^4  | (...)
            %                |------------------------------------------------------------------------------------------------|
            %  > Coeff     = | 1 | 1 | 1 | 1   | 1   | 1   | 1    | 1    | 1    | 1    | 1    | 1     | 1      | 1     | 1    | (...)
            %  > X_exp     = | 0 | 1 | 0 | 2   | 1   | 0   | 3    | 2    | 1    | 0    | 4    | 3     | 2      | 1     | 0    | (...)
            %  > Y_exp     = | 0 | 0 | 1 | 0   | 1   | 2   | 0    | 1    | 2    | 3    | 0    | 1     | 2      | 3     | 4    | (...)
            %                |------------------------------------------------------------------------------------------------|
            %                |------------------------------------------------------------------------------------------------|
            %  > dX        = | 0 | 1 | 0 | 2x  | y   | 0   | 3x^2 | 2xy  | y^2  | 0    | 4x^3 | 3x^2y | 2xy^2  | y^3   | 0    | (...)
            %  > Coeff     = | 0 | 1 | 0 | 2   | 1   | 0   | 3    | 2    | 1    | 0    | 4    | 3     | 2      | 1     | 0    | (...)
            %  > dX_exp(x) = | 0 | 0 | 0 | 1   | 0   | 0   | 2    | 1    | 0    | 0    | 3    | 2     | 1      | 0     | 0    | (...)
            %  > dX_exp(y) = | 0 | 0 | 0 | 0   | 1   | 0   | 0    | 1    | 2    | 0    | 0    | 1     | 2      | 3     | 0    | (...)
            %                |------------------------------------------------------------------------------------------------|
            %  > dY        = | 0 | 0 | 1 | 0   | x   | 2y  | 0    | x^2  | 2xy  | 3y^2 | 0    | x^3   | 2x^2y  | 3xy^2 | 4y^3 | (...)
            %  > Coeff     = | 0 | 0 | 1 | 0   | 1   | 2   | 0    | 1    | 2    | 3    | 0    | 1     | 2      | 3     | 4    | (...)
            %  > dX_exp(x) = | 0 | 0 | 0 | 0   | 1   | 0   | 0    | 2    | 1    | 0    | 0    | 3     | 2      | 1     | 0    | (...)
            %  > dX_exp(y) = | 0 | 0 | 0 | 0   | 0   | 1   | 0    | 0    | 1    | 2    | 0    | 0     | 1      | 2     | 3    | (...)
            %                |------------------------------------------------------------------------------------------------|
            %                |------------------------------------------------------------------------------------------------|
            %  > #         = | 1 | 2 | 3 | 4   | 5   | 6   | 7    | 8    | 9    | 10   | 11   | 12    | 13     | 14    | 15   | (...)
            %                |------------------------------------------------------------------------------------------------| (...)
            %                | p = 1     | p = 3                                | p = 5                                      ...
            %                |------------------------------------------------------------------------------------------------|
            
            % >> iD = 1.
            %  > Polynomial coefficients.
            Coeff_1 = ones(1,len_p);
            %  > Polynomial exponents.
            Exp_1   = zeros(2,len_p);
            i       = 1;
            while i < p+1
                %  > Previous last index.
                k = 1./2.*(i.^2+i);
                %  > X.
                for j = 1:i+1
                    Exp_1(1,k+j) = i+1-j;
                end
                %  > Y.
                for j = 1:i+1
                    Exp_1(2,k+j) = j-1;
                end
                i = i+1;
            end
            
            % >> iD = 2.
            %  > Polynomial exponents.
            if iD == 2
                [Exp_dX,Exp_dY] = deal(zeros(2,len_p));
                %  > grad_x(x): Polynomial exponents.
                i = 2;
                while i < p+1
                    %  > Previous last index.
                    k = 1./2.*(i.^2+i);
                    %  > X.
                    for j = 1:i
                        Exp_dX(1,k+j) = i-j;
                    end
                    %  > Y.
                    for j = 1:i
                        Exp_dX(2,k+j) = j-1;
                    end
                    i = i+1;
                end
                %  > grad_x(y): Polynomial exponents.
                i = 2;
                while i < p+1
                    %  > Previous last index.
                    k = 1./2.*(i.^2+i);
                    %  > X.
                    for j = 1:i
                        Exp_dY(1,k+j+1) = i-j;
                    end
                    %  > Y.
                    for j = 1:i
                        Exp_dY(2,k+j+1) = j-1;
                    end
                    i = i+1;
                end
            end
            
            % >> Deal coefficients/exponents.
            if iD == 1
                % >> Phi_f.
                Coeff_iD = Coeff_1;
                Exp_iD   = Exp_1;
            elseif iD == 2
                % >> grad(Phi_f).
                %  > grad(x).
                Coeff_iD (1,:) = Exp_1 (1,:);
                Exp_iD{1}(1,:) = Exp_dX(1,:);
                Exp_iD{1}(2,:) = Exp_dX(2,:);
                %  > grad(y).
                Coeff_iD (2,:) = Exp_1 (2,:);
                Exp_iD{2}(1,:) = Exp_dY(1,:);
                Exp_iD{2}(2,:) = Exp_dY(2,:);
            end
        end
        
        %% > 2. -----------------------------------------------------------
        function [wf_1,wf_2] = W_Function()
            % >> Symbolic variable.
            syms a b d;
            
            %  > Weight function(s).
            wf_1 = @(a,b,d) 1./(a.*d.^b);
            wf_2 = @(a,b,d) 1./(a.*(d./2).^b);
        end
    end
end