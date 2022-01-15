classdef SubClass_2_2
    methods (Static)
        %% > Wrap up SubClass_2_2.
        function [msh] = WrapUp_2_2(inp,msh,fn)
            % >> ----------------------------------------------------------
            % >> 1.   Assemble matrices Df, Dwf, Pf and Tf.
            %  > 1.1. Wrap up '2.' and call '2.2.'.
            %  > 1.2. Compute polynomial terms.
            %  > 1.3. Set weighting function.
            %  > 1.4. Compute df array.
            % >> 2.   Assemble matrices A and B.
            %  > 2.1. Assemble matrix A.
            %  > 2.2. Assemble matrix B.
            % >> ----------------------------------------------------------
            % >> Local variables.
            wf    = inp.fr.wf;
            ng    = inp.fr.ng;
            V     = [inp.pr.vx,inp.pr.vy];
            G     = [inp.pr.gx,inp.pr.gy];
                        
            [msh,Tf_C,Tf_D] = SubClass_2_2.WrapUp_2(msh,msh.s.xy_v,wf,ng);
            [msh] = SubClass_2_2.WrapUp_3(msh,V,G,Tf_C,Tf_D);
        end       
                
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [msh,Tf_C,Tf_D] = WrapUp_2(msh,st_xy,wf,ng)
            %% > Polynomial reconstruction.
            %  > Polynomial order.
            p = 2.*size(msh.s.st,1)-1;
            %  > Number of terms.
            len_p = 1./2.*(p+1).*(p+2);
            %  > Polynomial terms.
            [Coeff_1,Exp_1] = SubClass_2_2.Compute_PolTerms(1,p,len_p);
            [Coeff_2,Exp_2] = SubClass_2_2.Compute_PolTerms(2,p,len_p);
                 
            %% > Matrices Df and Dwf.
            % >> Df.
            %  > Df = [1*{(x-xf)^0}*{(y-yf)^0},1*{(x-xf)^1}*{(y-yf)^0},1*{(x-xf)^0}*{(y-yf)^1},...] = [1,(x-xf),(y-yf),...].
            for i = 1:msh.f.NF
                Face(i,:) = ...
                    [msh.f.mean(1,i),msh.f.mean(2,i)];
                for j = 1:size(st_xy{i},2)
                    XY{i}(1,j)  = st_xy{i}(1,j)-Face(i,1);
                    XY{i}(2,j)  = st_xy{i}(2,j)-Face(i,2);
                    for k = 1:len_p
                        Df{i}(j,k) = Coeff_1(k).*(XY{i}(1,j).^Exp_1(1,k)).*(XY{i}(2,j).^Exp_1(2,k));
                    end
                end
            end
            
            % >> Dwf.
            %  > Dwf = w*Df.
            %  > Boundary faces' cell index.
            for i = 1:size(msh.bnd.f,2)
                bnd_ff(i) = msh.bnd.f{2,i};
                bnd_fc(i) = msh.bnd.f{3,i};
            end
            for i = 1:msh.f.NF
                for j = 1:size(st_xy{i},2)
                    %  > Point.
                    %  > Remark: If point=face, then face "i" belongs to a corner cell -> use cell centroid instead.
                    if ~isequal(st_xy{i}(:,j)',Face(i,:))
                        Point{i}(j,:) = st_xy{i}(:,j)';
                        w    {i}(j,1) = SubClass_2_2.W_Function(Point{i}(j,:),Face(i,:),wf,4,'F');
                    else
                        Point{i}(j,:) = msh.c.mean(:,bnd_fc(bnd_ff == i));
                        w    {i}(j,1) = SubClass_2_2.W_Function(Point{i}(j,:),Face(i,:),wf,4,'T');
                    end
                    Dwf{i}(j,:) = Df{i}(j,:).*w{i}(j,1);
                end
            end
            
            %% > Matrices Pf and Tf. 
            % >> Pf.
            %  > Pf = inv(Dwf_T*Df)*Dwf_T.
            for i = 1:msh.f.NF
                Inv          {i} = eMatrices(transpose(Dwf{i})*Df{i});
                Pf           {i} = Inv{i}*transpose(Dwf{i});
                msh.s.cond_Df(i) = cond(Df{i});
            end           
            % >> Tf     = ?
            %  > C(1)   = P(1,1)*Phi(1)+P(1,2)*Phi(2)+P(1,3)*Phi(3)+P(1,4)*Phi(4)+...
            %  > C(2)   = P(2,1)*Phi(1)+P(2,2)*Phi(2)+P(2,3)*Phi(3)+P(2,4)*Phi(4)+...
            %  > C(3)   = P(3,1)*Phi(1)+P(3,2)*Phi(2)+P(3,3)*Phi(3)+P(3,4)*Phi(4)+...
            %  > Phi_f  = C(1)+C(2)*(x-xf)+C(3)*(y-yf)+...
            %  > Phi_f  = {P(1,1)*Phi(1)+P(1,2)*Phi(2)+P(1,3)*Phi(3)+P(1,4)*Phi(4)}+
            %           + {P(2,1)*Phi(1)+P(2,2)*Phi(2)+P(2,3)*Phi(3)+P(2,4)*Phi(4)}*(x-xf)+
            %           + {P(3,1)*Phi(1)+P(3,2)*Phi(2)+P(3,3)*Phi(3)+P(3,4)*Phi(4)}*(y-yf)+...
            %  > Phi_f  = {P(1,1)+P(2,1)*(x-xf)+P(3,1)*(y-yf)+...}*Phi(1)+
            %           + {P(1,2)+P(2,2)*(x-xf)+P(3,2)*(y-yf)+...}*Phi(2)+
            %           + {P(1,3)+P(2,3)*(x-xf)+P(3,3)*(y-yf)+...}*Phi(3)+...
            %           + {P(1,4)+P(2,4)*(x-xf)+P(3,4)*(y-yf)+...}*Phi(4)+...
            %
            %  > Tf    -> Equivalent to: [1,(x-xf),(y-yf),...]*Pf*Phi = df*Pf*Phi.
            %    e.g.:    p=1 w/ ns = 8 -> (1x3)*(3x8)*(8x1) = (1x8)*(8x1) = (1x1).
            
            % >> Face quadrature.
            %  > Function handle coordinate transformation.
            [xy_fg,j_fg] = SubClass_2_2.Create_FunctionHandle();
            for i = 1:msh.f.NF
                msh.f.gq{i} = SubClass_2_2.GaussFace_Points(xy_fg,j_fg,ng,msh.f.xy_v{i});
            end
            
            %  > Convection.
            for i = 1:msh.f.NF
                df  {i} = SubClass_2_2.Compute_df(1,Face(i,:),msh.f.gq{i},Coeff_1,Exp_1);
                Tf_C{i} = df{i}*Pf{i};
            end
            %  > Diffusion.
            for i = 1:msh.f.NF
                grad_df{i} = SubClass_2_2.Compute_df(2,Face(i,:),msh.f.gq{i},Coeff_2,Exp_2);
                Tf_D   {i} = grad_df{i}*Pf{i};
            end
        end
        % >> 2.2. ---------------------------------------------------------
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
        % >> 2.3. ---------------------------------------------------------
        function [W_func] = W_Function(Point,Face,wf,k,Flag)
            % >> Weighting functions: 1) Unweighted: wij = 1.
            %                         2) Weighted  : wij = 1/|xj-xi|^2 -> wij = 1/|dist(x_point-xf)|^k.
            %                         3) Weighted  : wij = (...).
            if strcmpi(wf,'Unweighted')
                %  > 1)
                W_func = 1;
            elseif strcmpi(wf,'Weighted')
                %  > 2)
                d = pdist([Point;Face],'euclidean');
                if strcmpi(Flag,'F')
                    W_func = 1./d.^k;
                elseif strcmpi(Flag,'T')
                    W_func = 1./(d./2).^k;
                    A = 1-exp(-k.^2);
                    W_func = 1;
                end
            end
        end
        % >> 2.4. ---------------------------------------------------------
        function [df_ij] = Compute_df(iD,mean_f,fg,Coeff,Exp)
            if iD == 1
                % >> Phi_f.
                for i = 1:size(fg.Points,1)
                    for j = 1:length(Coeff)
                        df_ij(i,j) = Coeff(j).*fg.Weights(i).*((mean_f(1)-fg.Points(i,1)).^Exp(1,j)).*((mean_f(2)-fg.Points(i,2)).^Exp(2,j));
                    end
                end
                df_ij = sum(df_ij,1);
            elseif iD == 2
                % >> grad(Phi_f).
                for i = 1:size(Coeff,1)
                    %  > i=1: X.
                    %  > i=2: Y.
                    for j = 1:size(fg.Points,1)
                        for k = 1:size(Coeff,2)
                            df_ij{i}(j,k) = Coeff(i,k).*fg.Weights(j).*((mean_f(1)-fg.Points(j,1)).^Exp{i}(1,k)).*((mean_f(2)-fg.Points(j,2)).^Exp{i}(2,k));
                        end
                    end
                    df_ij{i} = sum(df_ij{i},1);
                end
                %  > Reshape cell array.
                df_ij = cell2mat(reshape(df_ij,[2,1]));
            end
        end                      
    end
end