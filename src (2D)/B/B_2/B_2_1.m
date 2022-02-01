classdef B_2_1
    methods (Static)
        %% > 1. -----------------------------------------------------------
        % >> 1.1. ---------------------------------------------------------
        function [len_p,Coeff_1,Coeff_2,Exp_1,Exp_2] = Polynomial_Reconstruction(np_x,np_y)
            %  > Number of (maximum) terms.
            p     = max(np_x,np_y);
            len_p = 1./2.*(p+1).*(p+2);
            %  > Polynomial coefficients/exponents.
            [Coeff_1,Exp_1] = B_2_1.Compute_PolTerms(1,p,len_p,np_x,np_y);
            [Coeff_2,Exp_2] = B_2_1.Compute_PolTerms(2,p,len_p,np_x,np_y);
        end
        % >> 1.2. ---------------------------------------------------------
        function [Coeff_iD,Exp_iD] = Compute_PolTerms(iD,p,len_p,np_x,np_y)
            %% > Isotropy (x,y).
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
            
            if np_x == np_y
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
            
            %% > Anisotropy.
            %  > X-anisotropy: remove terms in y...
            %  > Y-anisotropy: remove terms in x...           
            if np_x ~= np_y
                %% > X.
                %  > Save terms in the linear profile and exclude higher-order ones.
                n_lin = 1+size(Exp_1,1);
                iY    = find(Exp_1(2,n_lin+1:len_p) > Exp_1(1,n_lin+1:len_p))+n_lin;
                iX    = setdiff(1:len_p,iY);
                if iD == 1
                    % >> Phi_f.
                    Coeff_iD = Coeff_1(:,iX);
                    Exp_iD   = Exp_1  (:,iX);
                elseif iD == 2
                    % >> grad(Phi_f).
                    %  > grad(x).
                    Coeff_iD (1,:) = Exp_1 (1,iX);
                    Exp_iD{1}(1,:) = Exp_dX(1,iX);
                    Exp_iD{1}(2,:) = Exp_dX(2,iX);
                    %  > grad(y).
                    Coeff_iD (2,:) = Exp_1 (2,iX);
                    Exp_iD{2}(1,:) = Exp_dY(1,iX);
                    Exp_iD{2}(2,:) = Exp_dY(2,iX);
                end
            end
        end
        
        %% > 2. -----------------------------------------------------------
        % >> 2.1. ---------------------------------------------------------
        function [wf_1,wf_2] = W_Function()
            % >> Symbolic variable.
            syms a b d;
            
            %  > Weight function(s).
            wf_1 = @(a,b,d) 1./(a.*d.^b);
            wf_2 = @(a,b,d) 1./(a.*(d./2).^b);
        end
        % >> 2.2. ---------------------------------------------------------
        function [df_ij] = Compute_df(iD,mean_f,fg,numb,Coeff,Exp)
            if iD == 1
                % >> Phi_f.
                i          = 1:size(fg.Points,1);
                j          = 1:numb;
                df_ij(i,j) = Coeff(j).*fg.Weights(i).*((mean_f(1)-fg.Points(i,1)).^Exp(1,j)).*((mean_f(2)-fg.Points(i,2)).^Exp(2,j));
                df_ij      = sum(df_ij,1);
            elseif iD == 2
                % >> gradPhi_f.
                for i = 1:size(Coeff,1)
                    j              = 1:size(fg.Points,1);
                    k              = 1:numb;
                    df_ijk{i}(j,k) = Coeff(i,k).*fg.Weights(j).*((mean_f(1)-fg.Points(j,1)).^Exp{i}(1,k)).*((mean_f(2)-fg.Points(j,2)).^Exp{i}(2,k));
                    df_ijk{i}      = sum(df_ijk{i},1);
                    df_ij (i,:)    = [df_ijk{i}];
                end
            end
        end
        % >> 2.3. ---------------------------------------------------------
        function [Tf_C,Tf_D] = Assemble_Tf(msh,np_x,np_y,ng,wf,a,b,bnd_fc,bnd_ff)
            %% > Df.
            % >> (x,y) coordinates of...
            i         = 1:msh.f.NF;
            %  > ...face centroid.
            Face(1,i) = msh.f.mean(1,i);
            Face(2,i) = msh.f.mean(2,i);
            %  > ...stencil points.
            stl_xy    = msh.s.xy_v_t;
            % >> Polynomial regression coefficients.
            [len_p,Coeff_1,Coeff_2,Exp_1,Exp_2] = B_2_1.Polynomial_Reconstruction(np_x,np_y);
            
            % >> Df = [1*{(x-xf)^0}*{(y-yf)^0},1*{(x-xf)^1}*{(y-yf)^0},1*{(x-xf)^0}*{(y-yf)^1},...] = [1,(x-xf),(y-yf),...].
            k     = 1:len_p;
            for i = 1:msh.f.NF
                %  > j.
                j {i}         = 1:size(stl_xy{i},2);
                %  > x-xf.
                XY{i}(1,j{i}) = stl_xy{i}(1,j{i})-Face(1,i);
                XY{i}(2,j{i}) = stl_xy{i}(2,j{i})-Face(2,i);
                %  > Polynomial coefficients.
                C1{i}   (:,k) = repelem(Coeff_1(k),length(j{i}),1);
                C2{i}   (:,k) = XY{i}(1,j{i})'.^Exp_1(1,k).*XY{i}(2,j{i})'.^Exp_1(2,k);
                Df{i}(j{i},k) = C1{i}.*C2{i};
                %  > Slower alternative:
                %  for k = 1:len_p
                %      Df{i}(j{i},k) = Coeff_1(k).*(XY{i}(1,j{i}).^Exp_1(1,k)).*(XY{i}(2,j{i}).^Exp_1(2,k));
                %  end
            end
            
            %% > Dwf.
            % >> Dwf = W*Df.
            if strcmpi(wf,'Unweighted')
                Dwf = Df;
            elseif strcmpi(wf,'Weighted')
                %  > Set weighting function.
                [wf_1,wf_2] = B_2_1.W_Function();
                
                for i = 1:msh.f.NF
                    %  > if d~=0...
                    d{i} = sqrt((stl_xy{i}(1,:)-Face(1,i)).^2+(stl_xy{i}(2,:)-Face(2,i)).^2);
                    %  > if d=0...
                    %    Remark: If point=face, use cell centroid instead.
                    d_flag{i} = find(~d{i});
                    if ~isempty(d_flag{i})
                        kf       = d_flag{i};
                        Point    = msh.c.mean(:,bnd_fc(bnd_ff == i));
                        d{i}(kf) = A_Tools.fft_dist([Point';Face(:,i)']);
                    end
                    %  > Check whether d=0...
                    if isempty(d_flag{i})
                        w  {i}(:,1)      = wf_1(a,b,d{i});
                    else
                        len              = length(d{i});
                        idf{i}           = setdiff(1:len,kf);
                        w  {i}(idf{i},1) = wf_1(a,b,d{i}(idf{i}));
                        w  {i}(kf    ,1) = wf_2(a,b,d{i}(kf));
                    end
                    Dwf{i} = Df{i}.*w{i};
                end
            else
                return;
            end
            
            %% > Pf.
            % >> Pf = inv(Dwf_T*Df)*Dwf_T.
            for i = 1:msh.f.NF
                %  > Df*W*Df.
                DwfDf{i} = transpose(Dwf{i})*Df{i}; %  > Square matrix.
                %  > Select...
                %  > See prof. Duarte's PhD thesis (page 56).
                if size(DwfDf{i},1) <= 4
                    Inv{i} = B_Tools.CramerRule(DwfDf{i});
                else
                    Inv {i} = DwfDf{i}\eye(size(DwfDf{i}));
                    %  > Alternatives:
                    %  [U{i},S{i},V{i}] = svd_lapack(DwfDf{i});
                    %  Inv{i}           = V{i}*(S{i}\U{i}');
                    %  Inv{i}           = inv(DwfDf{i});
                end
                Pf{i} = Inv{i}*transpose(Dwf{i});
            end
            
            %% > Tf.
            % >> Tf = [1,(x-xf),(y-yf),...]*Pf*Phi = df*Pf*Phi.
            %  > 1D Quadrature.
            [xy_fg,j_fg,Q_1D] = B_1_2.CD_1D(ng);
            
            for i = 1:msh.f.NF
                %  > Face quadrature.
                gq     {i} = B_1_2.GaussFace_Points(Q_1D,xy_fg,j_fg,msh.f.xy_v{i});
                %  > Convection/diffusion contribution(s).
                df     {i} = B_2_1.Compute_df(1,Face(:,i),gq{i},len_p,Coeff_1,Exp_1);
                grad_df{i} = B_2_1.Compute_df(2,Face(:,i),gq{i},len_p,Coeff_2,Exp_2);
                Tf_C   {i} = df{i}*Pf{i};
                Tf_D   {i} = grad_df{i}*Pf{i};
            end
        end
    end
end