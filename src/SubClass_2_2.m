classdef SubClass_2_2
    methods (Static)
        %% > Wrap up SubClass_2_2.
        function [msh] = WrapUp_2_2(inp,msh,fn)
            %  ------------------------------------------------------------
            % >> 1.   Set stencil coordinates (gather face/cell centroids' coordinates).
            %  > 1.1. Wrap up '1.' and call '1.2.'
            %  > 1.2. Find neighbouring faces' centroid coordinates.
            % >> 2.   Assemble matrices Df, Dwf, Pf and Tf.
            %  > 2.1. Wrap up '2.' and call '2.2.'.
            %  > 2.2. Compute polynomial terms.
            %  > 2.3. Set weighting function.
            %  > 2.4. Invert ill-conditioned matrix (from Filipe Diogo's thesis).
            %  > 2.5. Quadrature abcissas/weights.
            %  > 2.6. Coordinate transformation: csi->x/y.
            %  > 2.7. Compute df array.
            % >> 3.   Assemble matrices A and B.
            %  > 3.1. Assemble matrix A.
            %  > 3.2. Assemble matrix B.
            %
            %  > ----------------------------------------------------------
            % >> Local variables.
            wf = inp.fr.wf;
            ng = inp.fr.ng;
            V  = [inp.pr.vx,inp.pr.vy];
            G  = [inp.pr.gx,inp.pr.gy];
            
            % >> 1.
            msh = SubClass_2_2.WrapUp_1(msh);
            % >> 2.
            [msh,Tf_C,Tf_D] = SubClass_2_2.WrapUp_2(msh,wf,ng);
            % >> 3.
            %[msh] = SubClass_2_2.WrapUp_3(V,G,Tf_C,Tf_D);
            
        end
        
        %% > Tools.
        %% > 1.) ----------------------------------------------------------
        % >> 1.1.) --------------------------------------------------------
        function [msh] = WrapUp_1(msh)
            % >> iD's.
            %  > Boundary cells' index.
            for i = 1:size(msh.bnd.c,2)
                bnd_cc(i) = msh.bnd.c{2,i};
            end
            %  > Boundary faces' cell index.
            for i = 1:size(msh.bnd.f,2)
                bnd_ff(i) = msh.bnd.f{2,i};
                bnd_fc(i) = msh.bnd.f{3,i};
            end
            
            % >> Deal (Xv,Yv) coordinates of stencil elements.
            for i = 1:msh.f.NF
                for is = 1:size(msh.s.st,1)
                    %  > Stencil cell centroids.
                    len_c{i}(is) = length(msh.s.st{is,i});
                    for j = 1:len_c{i}(is)
                        st_v{is,i}(1,j) = msh.c.mean(1,msh.s.st{is,i}(j));
                        st_v{is,i}(2,j) = msh.c.mean(2,msh.s.st{is,i}(j));
                    end
                    %  > Check whether face i's stencil cells belong to the boundary. If so, add the respective face to the stencil.
                    Flag{is,i} = zeros(1,len_c{i}(is));
                    Flag{is,i} = ismembc(msh.s.st{is,i},bnd_cc);
                    if any(Flag{is,i} == 1)
                        %  > Call 1.2.)
                        st_v{is,i} = SubClass_2_2.Deal_FaceCoord(Flag{is,i},st_v{is,i},bnd_ff,bnd_fc,len_c{i}(is),msh.s.st{is,i},msh.f.mean);
                    end
                    msh.s.XY_v{is,i} = st_v{is,i};
                end
            end
        end
        % >> 1.2.) --------------------------------------------------------
        function [st_v] = Deal_FaceCoord(Flag,st_v,bnd_ff,bnd_fc,len_c,st,mean_f)
            %  > Add respective boundary faces.
            %  > Remark: a given boundary cell may contain more than 1 boundary face (see 3rd row of msh.bnd.f)
            j = 1;
            for i = 1:length(Flag)
                if Flag(i)
                    cf{j} = find(bnd_fc == st(i));
                    j     = j+1;
                end
            end
            %  > Faces to be added to the stencil.
            f_add = cell2mat(cf);
            for j = 1:length(f_add)
                st_v(1,j+len_c) = mean_f(1,bnd_ff(f_add(j)));
                st_v(2,j+len_c) = mean_f(2,bnd_ff(f_add(j)));
            end
        end
        
        %% > 2.) ----------------------------------------------------------
        % >> 2.1.) --------------------------------------------------------
        function [msh,Tf_C,Tf_D] = WrapUp_2(msh,wf,ng)
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
            %  > Loop through faces.
            for i = 1:msh.f.NF
                % >> Face centroid coordinates.
                Face(i,:) = ...
                    [msh.f.mean(1,i),msh.f.mean(2,i)];
                %  > Loop through stencil elements.
                for j = 1:size(msh.s.XY_v{i},2)
                    % >> Df     = [1*{(x-xf)^0}*{(y-yf)^0},1*{(x-xf)^1}*{(y-yf)^0},1*{(x-xf)^0}*{(y-yf)^1},...] = [1,(x-xf),(y-yf),...]
                    XY{i}(1,j)  = msh.s.XY_v{i}(1,j)-Face(i,1);
                    XY{i}(2,j)  = msh.s.XY_v{i}(2,j)-Face(i,2);
                    %  > Loop through number of polynomial coefficients/exponents.
                    for k = 1:len_p
                        Df{i}(j,k) = Coeff_1(k).*(XY{i}(1,j).^Exp_1(1,k)).*(XY{i}(2,j).^Exp_1(2,k));
                    end
                end
            end
            
            % >> Dwf.
            %  > Boundary faces' cell index.
            for i = 1:size(msh.bnd.f,2)
                bnd_ff(i) = msh.bnd.f{2,i};
                bnd_fc(i) = msh.bnd.f{3,i};
            end
            %  > Loop through faces.
            for i = 1:msh.f.NF
                %  > Loop through stencil elements.
                for j = 1:size(msh.s.XY_v{i},2)
                    %  > Point.
                    %  > Remark: If point=face, then face "i" belongs to a corner cell -> use cell centroid instead.
                    if ~isequal(msh.s.XY_v{i}(:,j)',Face(i,:))
                        Point{i}(j,:) = msh.s.XY_v{i}(:,j)';
                    else
                        Point{i}(j,:) = msh.c.mean(:,bnd_fc(bnd_ff == i));
                    end
                    w  {i}(j,1) = SubClass_2_2.W_Function(Point{i}(j,:),Face(i,:),wf);
                    Dwf{i}(j,:) = Df{i}(j,:).*w{i}(j,1);
                end
            end
            
            %% > Matrices Pf and Tf.           
            for i = 1:msh.f.NF
                %  > Pf = inv(Dwf_T*Df)*Dwf_T.
                Pf{i} = (transpose(Dwf{i})*Df{i})\transpose(Dwf{i});
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
                msh.f.fg{i} = SubClass_2_2.GaussFace_Points(xy_fg,j_fg,ng,msh.f.XY_v{i});
            end
            
            %  > Convection.
            for i = 1:msh.f.NF
                df  {i} = SubClass_2_2.Compute_df(1,Face(i,:),msh.f.fg{i},Coeff_1,Exp_1);
                Tf_C{i} = df{i}*Pf{i};
            end
            %  > Diffusion.
            for i = 1:msh.f.NF
                grad_df{i} = SubClass_2_2.Compute_df(2,Face(i,:),msh.f.fg{i},Coeff_2,Exp_2);
                Tf_D   {i} = grad_df{i}*Pf{i};
            end
        end
        % >> 2.2.) --------------------------------------------------------
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
        % >> 2.3.) --------------------------------------------------------
        function [W_func] = W_Function(Point,Face,wf)
            % >> Weighting functions: 1) Unweighted: wij = 1.
            %                         2) Weighted  : wij = 1/|xj-xi|^2 -> wij = 1/|dist(x_point-xf)|^k.
            %                         3) Weighted  : wij = (...).
            if strcmpi(wf,'Unweighted')
                %  > 1)
                W_func = 1;
            elseif strcmpi(wf,'Weighted')
                %  > 2)
                k      = 2;
                W_func = 1./pdist([Point;Face],'euclidean').^k;
            end
        end
        % >> 2.4.) --------------------------------------------------------
        function [X] = Invert_Matrix(Df,Dwf,Tol,iterMax)
            %  > Divide each column of Df by its maximum absolute value.
            for j = 1:size(Df,2)
                Ff  (j) = max(abs(Df(:,j)));
                Df(:,j) = Df(:,j)./Ff(j);
            end
            Af = sparse(transpose(Dwf)*Df);
            %  > Preconditioning.
            setup = struct('type','ilutp','droptol',1e-1);
            [L,U] = ilu(Af,setup);
            %  > Use bicgstabl solver for every column...
            B = zeros(size(Af,1),1);
            for j = 1:size(Af,2)
                B (j,1)    = 1;
                [X(:,j),~] = bicgstabl(Af,B,Tol,iterMax,L,U,[]);
            end
            %  > Divide each row of X by Ff(i).
            for i = 1:size(X,1)
                X(i,:) = X(i,:)./Ff(i);
            end
        end
        % >> 2.5.) --------------------------------------------------------
        function [fg] = GaussFace_Points(xy_fg,j_fg,ng,xy_f)
            %  > Gauss-Legendre quadrature (available in [Tools - Numerical]/[QuadTools]).
            Q_1D = quadGaussLegendre(ng);
            
            for i = 1:length(Q_1D.Points)
                fg.Points(i,1) = xy_fg(xy_f(1,1),xy_f(2,1),Q_1D.Points(i));
                fg.Points(i,2) = xy_fg(xy_f(1,2),xy_f(2,2),Q_1D.Points(i));    
            end
            fg.Weights = 1./2.*Q_1D.Weights;
            fg.jac     = [j_fg(xy_f(1,1),xy_f(1,2)),j_fg(xy_f(2,1),xy_f(2,2))];
        end
        % >> 2.6.) --------------------------------------------------------
        function [x,j] = Create_FunctionHandle()
             % >> Symbolic variables.
            syms a b csi;

            %  > x(csi) = a*(1-csi)/2+b*(1+csi)/2.
            x = a.*(1-csi)./2+b.*(1+csi)./2;
            x = matlabFunction(x);
            %  > j(csi) = d(x)/d(csi) = (b-a)/2.
            j = (b-a)./2;
            j = matlabFunction(j);
        end
        % >> 2.7.) --------------------------------------------------------
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
                       
        %% > 3.) ----------------------------------------------------------
        % >> 3.1.) --------------------------------------------------------
        function [] = WrapUp_3(V,G,Tf_C,Tf_D)
            % >> Assemble matrices A and B.
            %  > A.
            SubClass_2_2.Assemble(V,G,Tf_C,Tf_D);
            %  > B.
            SubClass_2_2.Assemble(V,G,Tf_C,Tf_D);
            
            
            
            
            
        end
        % >> 3.2.) --------------------------------------------------------
        function [] = Assemble_A(V,G,Tf_C,Tf_D)
            
            
            
            
        end
        % >> 3.3.) --------------------------------------------------------
        function [] = Assemble_B(V,G,Tf_C,Tf_D)
        end
    end
end